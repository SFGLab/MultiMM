"""Short facts about 3D chromatin, shown under the banner (one random fact per run)."""
import random

FACTS = (
    "Stretched out, the DNA of one human cell is about 2 metres long, yet it fits in a nucleus roughly 10 micrometres across.",
    "A nucleosome is about 147 base pairs of DNA wrapped around a core of eight histone proteins.",
    "Nucleosomes repeat roughly every 200 base pairs, giving the 'beads on a string' look of the chromatin fibre.",
    "Inside cells, chromatin is mostly an irregular, flexible fibre; the neat 30 nm fibre seen in test tubes is rare in vivo.",
    "Linker histone H1 binds the DNA between nucleosomes and helps compact the fibre.",
    "Hi-C counts how often pairs of genomic regions touch, and gave the first genome-wide maps of 3D folding in 2009.",
    "Hi-C maps show chromatin splitting into A compartments (active, gene-rich) and B compartments (inactive, gene-poor).",
    "Rao and colleagues (2014) refined compartments into subcompartments (A1, A2, B1, B2, B3) with distinct chromatin marks.",
    "Topologically associating domains (TADs) are typically a few hundred kilobases to about a megabase long in mammals.",
    "Loop extrusion: the ring-shaped cohesin complex reels DNA through itself until it is stopped at CTCF sites.",
    "CTCF anchors loops best when its two binding sites face each other, the so-called convergent rule.",
    "Cohesin extrudes DNA at roughly one kilobase per second in laboratory experiments.",
    "Removing cohesin makes TADs and loops vanish within hours, yet most genes change their activity only modestly.",
    "A Hi-C map of a human cell line shows around 10,000 chromatin loops.",
    "Enhancers can switch on genes hundreds of kilobases away by looping next to them.",
    "In developing limbs, the HoxD gene cluster flips between two TADs as cells move from proximal to distal fates.",
    "Mutations that delete a TAD boundary can rewire enhancers to the wrong gene and cause limb malformations.",
    "Each chromosome occupies its own region of the nucleus, a territory, and neighbours intermingle only at the edges.",
    "Small, gene-rich chromosomes (such as chr 19) tend to sit near the nuclear centre; large, gene-poor ones near the edge.",
    "Lamina-associated domains (LADs) are gene-poor stretches of 0.1 to 10 megabases that touch the nuclear lamina.",
    "The nuclear lamina is a mesh of lamin proteins lining the inside of the nuclear envelope, and it tethers heterochromatin.",
    "The nucleolus assembles around ribosomal DNA repeats on chromosomes 13, 14, 15, 21 and 22.",
    "Heterochromatin is dense and mostly silent, often marked by H3K9me3 and bound by HP1 proteins.",
    "Polycomb proteins deposit H3K27me3 and cluster silenced developmental genes into compact foci.",
    "HP1 proteins can condense into liquid-like droplets, hinting that phase separation helps build heterochromatin.",
    "Highly active genes often sit close to nuclear speckles, bodies rich in splicing factors.",
    "Early-replicating DNA largely coincides with the A compartment, and late-replicating DNA with the B compartment.",
    "The contact probability between two loci falls roughly as 1/distance along the genome over a range of megabases.",
    "The fractal globule is a knot-free way to fold a long polymer, which is why a Hilbert curve is a good starting shape.",
    "Condensin turns chromosomes into dense arrays of loops when a cell prepares to divide.",
    "During mitosis, TADs and compartments largely disappear as chromosomes condense into rods.",
    "Topoisomerase II cuts and rejoins DNA so that long chromosomes can untangle without knotting.",
    "Single-cell imaging shows that TAD-like domains vary from cell to cell; a Hi-C map is an average over many cells.",
    "Chromatin is never still: a tagged locus wanders in its neighbourhood, moving slowly and sub-diffusively.",
    "The inactive X chromosome folds into two huge domains and condenses into a compact body at the nuclear edge.",
    "In rod cells of nocturnal mammals, heterochromatin sits in the middle and acts like a lens that guides light.",
    "In olfactory neurons, the single active odorant-receptor gene joins an enhancer hub assembled from several chromosomes.",
    "Even bacteria fold their DNA: the 4.6 million base pairs of E. coli, about 1.5 mm long, are packed into a small nucleoid.",
    "Polymer models need very few ingredients (loops, attraction between like regions) to reproduce much of a Hi-C map.",
    "Cell-type-specific 3D folding is one reason the same genome can run very different programs in a neuron and a liver cell.",
)


def random_fact() -> str:
    """One fact chosen at random (a different one on most runs)."""
    return random.SystemRandom().choice(FACTS)
