# Sanitize SAM alignment flags and mate fields for HTSJDK ValidationStringency.STRICT
# compatibility after upstream contig/unmapped filtering and samtools fixmate.
#
# Standard POSIX awk (no gawk extensions required; 12-bit flag tables are
# precomputed once in BEGIN{} for O(1) lookup per alignment record).
#
# Modes:
#   Default (AUDIT_ONLY=0): emits fixed SAM header + alignment records to stdout
#                           and prints AUDIT_SUMMARY metrics to stderr.
#   Audit   (AUDIT_ONLY=1): prints only AUDIT_SUMMARY to stderr and exits:
#                             0 -> clean (total > 0 and violations == 0)
#                             2 -> violations detected (violations > 0)
#                             3 -> zero alignment records observed (total == 0)

BEGIN {
    FS = "\t";
    OFS = "\t";
    AUDIT_ONLY = (AUDIT_ONLY + 0);

    for (v = 0; v < 4096; v++) {
        t = v;
        b1   = t % 2; t = int(t / 2);
        b2   = t % 2; t = int(t / 2);
        b4   = t % 2; t = int(t / 2);
        b8   = t % 2; t = int(t / 2);
        b16  = t % 2; t = int(t / 2);
        b32  = t % 2; t = int(t / 2);
        b64  = t % 2; t = int(t / 2);
        b128 = t % 2;
        upper = v - (b1 + 2*b2 + 4*b4 + 8*b8 + 16*b16 + 32*b32 + 64*b64 + 128*b128);

        is_paired[v]     = b1;
        is_proper[v]     = b2;
        is_unmap[v]      = b4;
        has_mate_unm[v]  = b8;

        no_proper[v]     = v - 2 * b2;
        # Legal bits when unpaired (0x1=0): 0x4 (unmap), 0x10 (rev), plus 0x100|0x200|0x400|0x800
        clean_unp[v]     = 4 * b4 + 16 * b16 + upper;
        bad_p_proper[v]  = (b1 && b2 && (b4 || b8));
        bad_p_r1r2[v]    = (b1 && ((b64 + b128) != 1));
        force_unmap[v]   = (v - 2 * b2) + (b4 ? 0 : 4);
    }
}

/^@/ {
    if ($1 == "@SQ") {
        for (i = 2; i <= NF; i++) {
            if (substr($i, 1, 3) == "SN:") {
                sq[substr($i, 4)] = 1;
                break;
            }
        }
    }
    if (!AUDIT_ONLY) print;
    next;
}

{
    total++;
    f = $2 + 0;
    dirty = 0;

    if (f < 0 || f >= 4096) {
        f = ((f % 4096) + 4096) % 4096;
        dirty = 1;
        fixed_range++;
    }

    # Unplaced reads (RNAME == "*") must have UNMAP (0x4) set, no PROPER_PAIR (0x2), and POS == 0
    if ($3 == "*") {
        if (!is_unmap[f]) {
            f = force_unmap[f];
            dirty = 1;
            fixed_unplaced_flag++;
        }
        if ($4 != "0") {
            if (!AUDIT_ONLY) $4 = "0";
            dirty = 1;
            fixed_unplaced_pos++;
        }
    }

    # Unmapped reads (0x4) must never have PROPER_PAIR (0x2)
    if (is_unmap[f] && is_proper[f]) {
        f = no_proper[f];
        dirty = 1;
        fixed_unmapped_proper++;
    }

    # Case A: Unpaired read (0x1 == 0)
    if (!is_paired[f]) {
        clean_f = clean_unp[f];
        has_mc = (index($0, "\tMC:Z:") > 0);
        if (clean_f != f || $7 != "*" || $8 != "0" || $9 != "0" || has_mc) {
            fixed_unpaired++;
            if (!AUDIT_ONLY) {
                if (has_mc) {
                    out = $1 "\t" clean_f "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t*\t0\t0\t" $10 "\t" $11;
                    for (i = 12; i <= NF; i++) {
                        if (substr($i, 1, 5) != "MC:Z:") {
                            out = out "\t" $i;
                        }
                    }
                    print out;
                    next;
                }
                $2 = clean_f;
                $7 = "*";
                $8 = "0";
                $9 = "0";
                print;
            }
            next;
        }
        if (dirty && !AUDIT_ONLY) {
            $2 = f;
            print;
            next;
        }
        if (!AUDIT_ONLY) print;
        next;
    }

    # Case B: Paired read (0x1 == 1)
    # 1. PROPER_PAIR (0x2) is invalid if either read (0x4) or mate (0x8) is unmapped
    if (bad_p_proper[f]) {
        f = no_proper[f];
        dirty = 1;
        fixed_paired_proper++;
    }

    # 2. Check R1 XOR R2 and mate reference consistency
    bad_r1r2 = bad_p_r1r2[f];
    bad_rnext_contig = ($7 != "*" && $7 != "=" && !($7 in sq));
    bad_mate_mapped = (!has_mate_unm[f] && ($7 == "*" || $8 == "0"));
    bad_mate_unmap_pos = ($7 == "*" && $8 != "0");

    if (bad_r1r2 || bad_rnext_contig || bad_mate_mapped || bad_mate_unmap_pos) {
        if (bad_r1r2) fixed_bad_r1r2++;
        if (bad_rnext_contig) fixed_bad_rnext++;
        if (bad_mate_mapped) fixed_bad_matepos++;
        if (bad_mate_unmap_pos) fixed_bad_mate_unmap_pos++;

        if (bad_r1r2 || bad_rnext_contig || bad_mate_mapped) {
            clean_f = clean_unp[f];
            if (!AUDIT_ONLY) {
                has_mc = (index($0, "\tMC:Z:") > 0);
                if (has_mc) {
                    out = $1 "\t" clean_f "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t*\t0\t0\t" $10 "\t" $11;
                    for (i = 12; i <= NF; i++) {
                        if (substr($i, 1, 5) != "MC:Z:") {
                            out = out "\t" $i;
                        }
                    }
                    print out;
                    next;
                }
                $2 = clean_f;
                $7 = "*";
                $8 = "0";
                $9 = "0";
                print;
            }
            next;
        }

        # Only RNEXT=="*" with non-zero PNEXT/TLEN under mate-unmapped (0x8)
        if (!AUDIT_ONLY) {
            $2 = f;
            $8 = "0";
            $9 = "0";
            print;
        }
        next;
    }

    if (dirty && !AUDIT_ONLY) {
        $2 = f;
        print;
        next;
    }
    if (!AUDIT_ONLY) print;
}

END {
    violations = fixed_range + fixed_unplaced_flag + fixed_unplaced_pos + \
                 fixed_unmapped_proper + fixed_unpaired + fixed_paired_proper + \
                 fixed_bad_r1r2 + fixed_bad_rnext + fixed_bad_matepos + fixed_bad_mate_unmap_pos;

    printf("AUDIT_SUMMARY\ttotal=%d\tviolations=%d\tfixed_unpaired=%d\tfixed_unmapped_proper=%d\tfixed_unplaced_pos=%d\tfixed_unplaced_flag=%d\tfixed_paired_proper=%d\tfixed_bad_r1r2=%d\tfixed_bad_rnext=%d\tfixed_bad_matepos=%d\tfixed_bad_mate_unmap_pos=%d\tfixed_range=%d\n",
        total + 0,
        violations + 0,
        fixed_unpaired + 0,
        fixed_unmapped_proper + 0,
        fixed_unplaced_pos + 0,
        fixed_unplaced_flag + 0,
        fixed_paired_proper + 0,
        fixed_bad_r1r2 + 0,
        fixed_bad_rnext + 0,
        fixed_bad_matepos + 0,
        fixed_bad_mate_unmap_pos + 0,
        fixed_range + 0) > "/dev/stderr";

    if (AUDIT_ONLY) {
        if ((total + 0) == 0) exit 3;
        if (violations > 0) exit 2;
    }
}
