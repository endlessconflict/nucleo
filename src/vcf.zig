//! VCF site records (CHROM, POS, REF, ALT) from a buffer already in memory,
//! with no allocation. Header lines are skipped; genotype columns are ignored.

const std = @import("std");

pub const Record = struct {
    chrom: []const u8,
    /// 1-based, as in the file.
    pos: u64,
    ref: []const u8,
    /// Comma-separated alternate alleles; iterate with `alts()`.
    alt: []const u8,

    pub fn alts(r: Record) std.mem.SplitIterator(u8, .scalar) {
        return std.mem.splitScalar(u8, r.alt, ',');
    }
};

pub const Error = error{BadRecord};

pub const Parser = struct {
    buf: []const u8,
    pos: usize = 0,

    pub fn next(self: *Parser) Error!?Record {
        while (self.pos < self.buf.len) {
            const end = std.mem.indexOfScalarPos(u8, self.buf, self.pos, '\n') orelse self.buf.len;
            var l = self.buf[self.pos..end];
            self.pos = @min(end + 1, self.buf.len);
            if (l.len > 0 and l[l.len - 1] == '\r') l = l[0 .. l.len - 1];
            if (l.len == 0 or l[0] == '#') continue;
            var f = std.mem.splitScalar(u8, l, '\t');
            const chrom = f.next() orelse return error.BadRecord;
            const pos = std.fmt.parseInt(u64, f.next() orelse return error.BadRecord, 10) catch return error.BadRecord;
            _ = f.next() orelse return error.BadRecord; // ID
            const ref = f.next() orelse return error.BadRecord;
            const alt = f.next() orelse return error.BadRecord;
            if (pos == 0 or ref.len == 0 or alt.len == 0) return error.BadRecord;
            return .{ .chrom = chrom, .pos = pos, .ref = ref, .alt = alt };
        }
        return null;
    }
};

test "sites, multi-allelic ALT, header and CRLF" {
    const text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\n" ++
        "chr22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\r\n" ++
        "chr22\t105\trs1\tAC\tA,ACC\t50\n";
    var p = Parser{ .buf = text };
    const a = (try p.next()).?;
    try std.testing.expectEqualStrings("chr22", a.chrom);
    try std.testing.expectEqual(@as(u64, 100), a.pos);
    try std.testing.expectEqualStrings("G", a.alt);
    const b = (try p.next()).?;
    var it = b.alts();
    try std.testing.expectEqualStrings("A", it.next().?);
    try std.testing.expectEqualStrings("ACC", it.next().?);
    try std.testing.expect(it.next() == null);
    try std.testing.expect(try p.next() == null);
}

test "malformed line is an error" {
    var p = Parser{ .buf = "chr1\tx\t.\tA\tC\n" };
    try std.testing.expectError(error.BadRecord, p.next());
}
