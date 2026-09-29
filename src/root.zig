//! nucleo: sequence I/O and k-mer primitives for the Zig bio libraries.

pub const fastx = @import("fastx.zig");
pub const kmer = @import("kmer.zig");
pub const vcf = @import("vcf.zig");

test {
    _ = fastx;
    _ = kmer;
    _ = vcf;
}
