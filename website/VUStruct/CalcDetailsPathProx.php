<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>PathProx</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">
    <link href="./css/vustruct.css" rel="stylesheet">
</head>
<body>
<div class="container">
    <?php echo file_get_contents("./navbar.html") ?>
<h2>PathProx: Pathogenic Proximity</h2>
 <p>Ranges of PathProx scores are presented for each variant on the generated case-wide home page.
    The individual ΔΔG calculations for each structure are in the Structure Summary near the top
    of each variant-specific drill-down page.  We now describe how we interpret these values.</p>
<p>
Developed by former Ph.D. student Mike Sivley in the Capra Lab, PathProx(<a
        href="https://doi.org/10.1016/j.ajhg.2018.01.017" target="_blank">Sivley et al. 2018</a>, <a
        href="https://doi.org/10.1186/s12859-018-2010-z" target="_blank">Sivley et al. 2018</a>)
    analyzes the placement of Clinvar(<a
        href="https://doi.org/10.1093/nar/gkx1153" target="_isblank">Landrum et al. 2018</a>),
    Gnomad(<a
        href="https://doi.org/10.1038/s41586-020-2308-7" target="_blank">Karczewski et al. 2020</a>),
    and random (null hypothesis) genetic variants on pipeline-selected 3D structures, and calculates whether the
    location of patient’s VUS better clusters with the pathogenic or benign variant sets.
    (Mike also authored the precursor to today’s VUStruct pipeline).
</p>
<p>
The AJHG citation(<a
        href="https://doi.org/10.1016/j.ajhg.2018.01.017" target="_blank">Sivley et al. 2018</a>) offers a thorough
    explanation of the algorithm’s mathematics and machine-learning approach.  The intuitive rationale for the
    algorithm comes from the observation that, often, when we map variants onto amino acid side chains,
    we observe a higher percentage of “known pathogenic” variants clustering inside smaller spheres vs.
    randomly distributed variants.  This technique identifies 3D hot-spots in the cartesian space of the folded protein.
</p>
<p>
For each transcript PathProx mines the Clinvar database for variants annotated as “likely pathogenic” or “pathogenic.”
    Benign variants are harvested from the Genome Aggregation Database gnomAD and filtered for MAF > 1x10<sup>-5</sup>.
    Roughly, we only accept as benign those variants which are found in at least 2 gnomAD sequences,
    a threshold that well contrasts to the rarity of UDN patient variants.
</p>
<p>
As a control, the PathProx algorithm computes a null distribution via thousands of random distributions of both
N=(count Gnomad) benign and M=(count Clinvar) pathogenic spheres on the then asks:
“Relative to random positioning of N spheres, how do these specific N pathogenic (or benign) variants
seem to cluster?”.  Often, we see significantly tighter clustering of the known pathogenic variants than
would be expected from our random sampling.  Simultaneously we see more diffuse clustering of benign variants vs.
random placements.  When we have both observations, and when our VUS fits more closely with the tighter
clustering pathogenic set, we identify it as pathogenic, and look more closely at the 3D context of the variant sets.
</p>
<p>
PathProx cannot always make a prediction and for these cases we leave a blank cell in our summary spreadsheet.
    PathProx minimally requires 3 pathogenic variants, and often leave-one-out validation tells us (via low AUC),
    that the predictive power of the algorithm is not strong.  When the calculation metrics imply confidence,
    we add a “no”, “maybe”, or “yes” to the PathProx-Clinvar column of the final spreadsheet.
</p>
<p>
As with the ΔΔG calculations, the PathProx algorithm is perhaps best validated across large protein datasets(<a
    href="https://doi.org/10.1016/j.ajhg.2018.01.017" target="_blank">Sivley et al. 2018</a>)
    Nonetheless, a strong case has been made for PathProx's power to explain disease in individual patients(<a
        href="https://doi.org/10.1186/s12859-018-2010-z" target="_blank">Sivley et al. 2018</a>)
</p>
<p>
As with ΔΔG, the PathProx algorithm was benchmarked on protein monomers.
    We nonetheless run it also on multimers from Swiss-model and the PDB.
    We cautiously consider these multimeric calculations in our analysis – as homomeric structures introduce symmetries
    of variant placements that the algorithm may not be adept at scoring.
</p>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4" crossorigin="anonymous"></script>
<script src="js/includeHTML.js">
</script>

<script>
    includeHTML();
</script>
</body>
</html>