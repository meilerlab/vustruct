<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Evolutionary Constraint</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">
    <link href="./css/vustruct.css" rel="stylesheet">
</head>
<body>
<div class="container">
    <?php echo file_get_contents("./navbar.html") ?>
<h2>Evolutionary Constraint: CONSURF and COSMIS</h2>
<p>
Residue positions which vary infrequently across species exhibit evolutionary conservation and this is already an input
to the clinic via GERP(<a
    href="https://doi.org/10.1101/gr.3577405" target="_blank">Cooper et al. 2005</a>) and GERP++(<a
    href="https://doi.org/10.1371/journal.pcbi.1001025" target="_blank">Davydov et al. 2010</a>) scores.
</p><p>
As structural biologists, we are interested not only in conservation for specific amino acids in a sequence, but also
in the conservation, and visualization, of the surrounding 3D protein context. ConSurf(<a
    href="https://doi.org/10.1093/nar/gkw408" target="_blank">Ashkenazy et al. 2016</a>) colors every
residue based on rates of evolution from phylogenetic tree analysis (<a
    href="https://doi.org/10.1093/bioinformatics/18.suppl_1.S71" target="_blank">Pupko et al. 2002</a> ).
</p>
<p>
The Contact Set MISsense tolerance algorithm, COSMIS (<a
        href="https://doi.org/10.1038/s41467-022-30936-x" target="_blank">Li et al. 2022</a>) was developed by Bian Li
of the Capra Lab. In contrast to Consurf, which quantifies evolutionary constraint based on analysis of sequence
alignments across diverse species, COSMIS quantifies evolutionary constraint over more recent evolution by
analyzing patterns of genetic variation across diverse human populations. It also leverages the 3D spatial context of
proteins to better estimate the evolutionary constraint on protein positions of interest. The figure below is an
example depiction of COSMIS scores.
</p>
    <figure class="center">
    <img src="./COSMIScoloringScreenShot.png"/>
        <figcaption>COSMIS scores are depicted on the protein backbone.
            Red communicates high conservation.
            Blue tolerance to variation.
            This coloring scheme is also available for scores from CONSURF and PathProx.  The colors in the example above reflect a general observation, that residues in highly structured protein regions tend to be more conserved across sequences, whether compared intra- or inter-species.  Our pipeline additionally can rainbow-color chains in complex structures, and show Alphafold pLDDTs with standard colors for those confidence metrics.</figcaption>
    </figure>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4" crossorigin="anonymous"></script>
<script src="js/includeHTML.js">
</script>

<script>
    includeHTML();
</script>
</body>
</html>