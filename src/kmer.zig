//! 2-bit nucleotide codes and canonical k-mers.
//!
//! A=0, C=1, G=2, T=3, so the complement of a code x is 3 - x (= x ^ 3), and
//! the reverse complement of a packed k-mer is a bit reversal of 2-bit groups
//! followed by XOR with all-ones. Any other byte (N, IUPAC) is invalid and
//! restarts k-mer extraction.

const std = @import("std");

pub const invalid: u8 = 0xff;

pub const code: [256]u8 = blk: {
    var t = [_]u8{invalid} ** 256;
    for ("ACGT", 0..) |c, i| {
        t[c] = i;
        t[std.ascii.toLower(c)] = i;
    }
    break :blk t;
};

/// Reverse complement of a k-mer packed with the first base in the high bits.
pub fn revcomp(x: u64, k: u6) u64 {
    // Reverse the order of 2-bit groups, then complement every base.
    var y = x;
    y = (y >> 2 & 0x3333333333333333) | (y & 0x3333333333333333) << 2;
    y = (y >> 4 & 0x0f0f0f0f0f0f0f0f) | (y & 0x0f0f0f0f0f0f0f0f) << 4;
    y = @byteSwap(y);
    return ~y >> @intCast(64 - 2 * @as(u7, k));
}

pub const Kmer = struct { pos: usize, fwd: u64, canonical: u64 };

/// Every valid k-mer (1 <= k <= 32) of `seq`, forward and canonical, skipping
/// windows that contain a non-ACGT byte. O(1) per base: rolling update of
/// forward and reverse-complement codes.
pub const Iterator = struct {
    seq: []const u8,
    k: u6,
    i: usize = 0,
    valid: usize = 0,
    fwd: u64 = 0,
    rev: u64 = 0,

    pub fn init(seq: []const u8, k: u6) Iterator {
        std.debug.assert(k >= 1 and k <= 32);
        return .{ .seq = seq, .k = k };
    }

    pub fn next(self: *Iterator) ?Kmer {
        const shift: u7 = 2 * @as(u7, self.k);
        const mask: u64 = if (shift == 64) ~@as(u64, 0) else (@as(u64, 1) << @intCast(shift)) - 1;
        while (self.i < self.seq.len) {
            const c = code[self.seq[self.i]];
            self.i += 1;
            if (c == invalid) {
                self.valid = 0;
                continue;
            }
            self.fwd = (self.fwd << 2 | c) & mask;
            self.rev = self.rev >> 2 | @as(u64, 3 - c) << @intCast(shift - 2);
            self.valid += 1;
            if (self.valid >= self.k) return .{
                .pos = self.i - self.k,
                .fwd = self.fwd,
                .canonical = @min(self.fwd, self.rev),
            };
        }
        return null;
    }
};

fn pack(s: []const u8) u64 {
    var x: u64 = 0;
    for (s) |c| x = x << 2 | code[c];
    return x;
}

test "revcomp of packed k-mers" {
    try std.testing.expectEqual(pack("ACGTT"), revcomp(pack("AACGT"), 5));
    try std.testing.expectEqual(pack("GGGC"), revcomp(pack("GCCC"), 4));
    const k32 = "ACGTACGTTTGCAGGCATCGATCGGATCAGTC";
    var rc: [32]u8 = undefined;
    for (k32, 0..) |c, i| rc[31 - i] = "TGCA"[code[c]];
    try std.testing.expectEqual(pack(&rc), revcomp(pack(k32), 32));
}

test "iterator: rolling codes equal direct packing, N restarts the window" {
    const seq = "ACGTNNACGGTTCAnACGTAC";
    const k = 4;
    var it = Iterator.init(seq, k);
    var count: usize = 0;
    while (it.next()) |km| : (count += 1) {
        const w = seq[km.pos..][0..k];
        for (w) |c| try std.testing.expect(code[c] != invalid);
        try std.testing.expectEqual(pack(w), km.fwd);
        try std.testing.expectEqual(@min(pack(w), revcomp(pack(w), k)), km.canonical);
    }
    // Valid windows: "ACGT" (1) + "ACGGTTCA" (5) + "ACGTAC" (3).
    try std.testing.expectEqual(@as(usize, 9), count);
}
