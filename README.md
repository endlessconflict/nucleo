# nucleo

Sequence input and k-mer primitives for Zig, kept small on purpose. It is the shared base for a set of bioinformatics libraries in Zig, and it only grows when one of them needs something.

## What is in it

`fastx` parses FASTA and FASTQ from a buffer that is already in memory, with no allocation. The parser compacts multi-line FASTA records in place, squeezing the newlines out of a record by moving bytes left within that record's own span, so each sequence you get back is one contiguous slice. However, this rewrites the buffer, so parse a copy if you still need the original text. Windows line endings and blank lines between records are accepted. FASTQ is read in its four-line form, and a quality line whose length differs from the sequence is an error rather than a silent truncation.

`kmer` maps bases to 2-bit codes (A=0, C=1, G=2, T=3) and walks every k-mer of a sequence for k up to 32, returning both the forward code and the canonical one (the smaller of forward and reverse complement). Each base updates both codes in constant time. Any byte other than ACGT, such as N, restarts the window. The reverse complement of a packed k-mer is a reversal of 2-bit groups followed by a bitwise NOT, which the code does with two mask-and-shift steps and a byte swap.

## Usage

Requires Zig 0.16.0.

```zig
const nucleo = @import("nucleo");

var p = nucleo.fastx.FastaParser{ .buf = buffer }; // buffer: []u8
while (try p.next()) |rec| {
    var it = nucleo.kmer.Iterator.init(rec.seq, 21);
    while (it.next()) |km| {
        _ = km.canonical;
    }
}
```

```sh
zig build test
```

The tests include a randomized comparison of the FASTA parser against a naive line-splitting parser over 2000 generated files.

## Not here yet

Multi-line FASTQ, compressed input, alignment and variant formats, and minimizer schemes. Each one lands when a library that depends on nucleo needs it.

## License

MIT
