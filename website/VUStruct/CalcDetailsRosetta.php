<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Rosetta ΔΔG</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">
    <link href="./css/vustruct.css" rel="stylesheet">
</head>
<body>
<div class="container">
    <?php echo file_get_contents("./navbar.html") ?>
<h2>Rosetta ΔΔG<sub>folding</sub></h2>
<p>
    Ranges of ΔΔG estimates are presented for each gene on the generated case-wide home page.
    The individual ΔΔG calculations for each structure are in the Structure Summary near the top
    of each variant-specific page.  We now describe how we interpret these values.
</p><p>
Proteins are tolerant of (i.e., are functional after) most single random amino acid substitutions(<a
        href="https://doi.org/10.1073/pnas.0403255101" target="_blank">Guo et al. 2004<a>,
        <a href="https://doi.org/10.1126/science.2315699" target="_blank">Bowie et al. 1990</a>).
    However, protein folding is an extraordinarily subtle phenomenon of nature.
    For globular protein, ∆G<sub>unfolded&rarr;folded</sub> measurements have long been understood to lie in the range [-5, -15] kCal/mol(<a
                    href="https://doi.org/10.1016/0968-0004(90)90124-T" target="_blank">Pace,C.N. 1990</a>).
    Thus, the free energy of protein folding is comparable to the enthalpy of a couple of hydrogen bonds - remarkably small relative to the sum
    of energetics acting between all the atoms of a protein or the global loss of entropy on folding.
    From basic thermodynamics, it must follow that a 1.5 kCal/mol increase to ΔG<sub>folding</sub> drives a 10-fold reduction in
    folded protein concentration and calculations have suggested that
    20% of deleterious missense variants cross this low energetic barrier(<a
            href="https://doi.org/10.1371/journal.pone.0107353" target="_blank">Berliner et al. 2014</a>).
    Though the prevailing theoretical intuition around protein folding thermodynamics is arguably oversimplistic
    when considering more complex proteins and in-vivo chaperones (<a
            href="https://doi.org/10.3390/ijms23010521" target="_blank">Sorokina et al. 2022</a>), the energetic impact
    of a single amino acid variant can obviously be sufficient to disrupt folding.  Trivial mechanistic explanations
    might include disruption of a critical hydrogen bond, loss of a disulfide bride,
    a new steric clash – any of which has the potential to render a protein unstable and disrupt function.
    Importantly, our calculations do not claim to calculate the thermodynamics of ∆G.  Rather, our alchemistic hope
    is that the comparison of varying energetics before and after variation will be sufficient to estimate a change
    in ∆G, a Δ(ΔG<sub>unfolded&rarr;folded</sub>)<sub>variant</sub>.
</p><p>
In comparison to the subtlety of nature, our Rosetta ΔΔG Cartesian (<a
        href="https://doi.org/10.1021/acs.jctc.6b00819" target="_blank">Park et al. 2016</a>, <a
        href="https://doi.org/10.3389/fbioe.2020.558247" target="_blank">Frenz et al. 2020</a>) calculations are a fast approximation.
    From a background total computed energy which includes hundreds or thousands of amino acid atomic interactions,
    (easily summing to 100,000s of kcal/mol) the Rosetta calculation attempts to determine how the new variant side
    chain would fit in a still-folded protein, allowing for limited rearrangement of nearby amino acids.
    The approximation of general protein rigidity, the treatment of only protein monomers without bound ligands, and
    the lack of membrane context make for a convenient and fast-running calculation with a result that we do not
    trust blindly.  Indeed, detailed benchmarking of the computed ΔΔGs on “simple” soluble monomeric proteins
    found its performance only somewhat helpful in the task of engineering large libraries of protein designs, an
    application quite distanced from characterizing individual protein variants in individualized patients(<a
        href="https://doi.org/10.3389/fbioe.2020.558247" target="_blank">Frenz et al. 2020</a>).
</p><p>
We nonetheless find the ΔΔG calculation illuminating in a reverse sense.  When variants have a small calculated energetic impact, in the interval [-2,2] kCal/mol for or purposes, we will generally classify these as benign, depending on structural context.  Outside that small range, we look more closely at 3D structures interactively, and we admit to layering subjectivity on the interpretation and classification of these results.  We discount the predictive power of the calculation when it involves residues from low-confidence models, or poorly resolved experimental structure.  Often, we find that the ΔΔG calculation will simply highlight amino acid changes that are predictably deleterious without any detailed calculation (ex: to or from proline or glycine, introduction of tryptophan, changes from large to small, hydrophobicity, etc.).  In the end, red/yellow/green are assigned to the spreadsheet cell after some reflection.  More robust analysis techniques are available if there is interest from the clinical team.
</p>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4" crossorigin="anonymous"></script>
<script src="js/includeHTML.js">
</script>

<script>
    includeHTML();
</script>
</body>
</html>