//! FASTA and FASTQ parsing over an in-memory buffer, with no allocation.
//!
//! Multi-line FASTA sequences are compacted in place: the newlines inside a
//! record are squeezed out by moving bytes left within that record's own span,
//! so every returned sequence is one contiguous slice of the input buffer.
//! Parsing therefore rewrites the buffer. Parse a copy if you need the original.

const std = @import("std");

pub const Record = struct {
    name: []const u8,
    seq: []const u8,
    /// FASTQ only.
    qual: ?[]const u8 = null,
};

pub const Error = error{ BadHeader, TruncatedRecord, QualityLengthMismatch };

/// Line starting at `pos` without its terminator (`\n` or `\r\n`), and the
/// position just past the terminator.
fn line(buf: []const u8, pos: usize) struct { []const u8, usize } {
    const end = std.mem.indexOfScalarPos(u8, buf, pos, '\n') orelse buf.len;
    const next = @min(end + 1, buf.len);
    const trimmed = if (end > pos and buf[end - 1] == '\r') end - 1 else end;
    return .{ buf[pos..trimmed], next };
}

pub const FastaParser = struct {
    buf: []u8,
    pos: usize = 0,

    pub fn next(self: *FastaParser) Error!?Record {
        // Skip blank lines before a record.
        while (self.pos < self.buf.len and (self.buf[self.pos] == '\n' or self.buf[self.pos] == '\r')) self.pos += 1;
        if (self.pos >= self.buf.len) return null;
        if (self.buf[self.pos] != '>') return error.BadHeader;
        const header, var pos = line(self.buf, self.pos + 1);
        const start = pos;
        var len: usize = 0;
        while (pos < self.buf.len and self.buf[pos] != '>') {
            const chunk, pos = line(self.buf, pos);
            std.mem.copyForwards(u8, self.buf[start + len ..][0..chunk.len], chunk);
            len += chunk.len;
        }
        self.pos = pos;
        return .{ .name = header, .seq = self.buf[start..][0..len] };
    }
};

/// Four-line FASTQ. ponytail: multi-line FASTQ (legal but rare) is rejected
/// as QualityLengthMismatch; support it when real data needs it.
pub const FastqParser = struct {
    buf: []const u8,
    pos: usize = 0,

    pub fn next(self: *FastqParser) Error!?Record {
        while (self.pos < self.buf.len and (self.buf[self.pos] == '\n' or self.buf[self.pos] == '\r')) self.pos += 1;
        if (self.pos >= self.buf.len) return null;
        if (self.buf[self.pos] != '@') return error.BadHeader;
        const header, var pos = line(self.buf, self.pos + 1);
        if (pos >= self.buf.len) return error.TruncatedRecord;
        const seq, pos = line(self.buf, pos);
        if (pos >= self.buf.len or self.buf[pos] != '+') return error.TruncatedRecord;
        _, pos = line(self.buf, pos);
        if (pos >= self.buf.len and seq.len > 0) return error.TruncatedRecord;
        const qual, pos = line(self.buf, pos);
        if (qual.len != seq.len) return error.QualityLengthMismatch;
        self.pos = pos;
        return .{ .name = header, .seq = seq, .qual = qual };
    }
};

// ---------------------------------------------------------------- tests

/// Obvious line-splitting reference parser for differential testing.
fn naiveFasta(gpa: std.mem.Allocator, text: []const u8, names: *std.ArrayList([]const u8), seqs: *std.ArrayList(std.ArrayList(u8))) !void {
    var it = std.mem.splitScalar(u8, text, '\n');
    while (it.next()) |raw| {
        const l = std.mem.trimEnd(u8, raw, "\r");
        if (l.len > 0 and l[0] == '>') {
            try names.append(gpa, l[1..]);
            try seqs.append(gpa, .empty);
        } else if (seqs.items.len > 0) {
            try seqs.items[seqs.items.len - 1].appendSlice(gpa, l);
        }
    }
}

test "FASTA: multi-line, CRLF, blank lines, empty sequence" {
    var text = ">a desc\r\nACGT\r\nGG\n\n>b\n>c\nNNNN\nA".*;
    var p = FastaParser{ .buf = &text };
    const a = (try p.next()).?;
    try std.testing.expectEqualStrings("a desc", a.name);
    try std.testing.expectEqualStrings("ACGTGG", a.seq);
    const b = (try p.next()).?;
    try std.testing.expectEqualStrings("b", b.name);
    try std.testing.expectEqualStrings("", b.seq);
    const c = (try p.next()).?;
    try std.testing.expectEqualStrings("NNNNA", c.seq);
    try std.testing.expect(try p.next() == null);
}

test "FASTA: randomized differential test against the naive parser" {
    const gpa = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(42);
    const rnd = prng.random();
    for (0..2000) |_| {
        var text: std.ArrayList(u8) = .empty;
        defer text.deinit(gpa);
        const nrec = rnd.intRangeAtMost(usize, 0, 6);
        for (0..nrec) |r| {
            try text.print(gpa, ">rec{d} x\n", .{r});
            for (0..rnd.intRangeAtMost(usize, 0, 5)) |_| {
                for (0..rnd.intRangeAtMost(usize, 0, 12)) |_| try text.append(gpa, "ACGTN"[rnd.intRangeLessThan(usize, 0, 5)]);
                try text.appendSlice(gpa, if (rnd.boolean()) "\r\n" else "\n");
            }
        }
        var names: std.ArrayList([]const u8) = .empty;
        defer names.deinit(gpa);
        var seqs: std.ArrayList(std.ArrayList(u8)) = .empty;
        defer {
            for (seqs.items) |*s| s.deinit(gpa);
            seqs.deinit(gpa);
        }
        try naiveFasta(gpa, text.items, &names, &seqs);

        const copy = try gpa.dupe(u8, text.items);
        defer gpa.free(copy);
        var p = FastaParser{ .buf = copy };
        var i: usize = 0;
        while (try p.next()) |rec| : (i += 1) {
            try std.testing.expectEqualStrings(names.items[i], rec.name);
            try std.testing.expectEqualStrings(seqs.items[i].items, rec.seq);
        }
        try std.testing.expectEqual(names.items.len, i);
    }
}

test "FASTQ: '@' inside quality, CRLF, errors" {
    var p = FastqParser{ .buf = "@r1\nACGT\n+\n@@II\r\n@r2\n\n+r2\n\n" };
    const r1 = (try p.next()).?;
    try std.testing.expectEqualStrings("ACGT", r1.seq);
    try std.testing.expectEqualStrings("@@II", r1.qual.?);
    const r2 = (try p.next()).?;
    try std.testing.expectEqualStrings("", r2.seq);
    try std.testing.expect(try p.next() == null);

    var bad = FastqParser{ .buf = "@r\nACGT\n+\nII\n" };
    try std.testing.expectError(error.QualityLengthMismatch, bad.next());
    var trunc = FastqParser{ .buf = "@r\nACGT\n" };
    try std.testing.expectError(error.TruncatedRecord, trunc.next());
}
