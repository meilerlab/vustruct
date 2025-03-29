<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Calculations Overview</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">
    <link href="./css/vustruct.css" rel="stylesheet">
</head>
<body>
<div class="container">
    <?php echo file_get_contents("./navbar.html") ?>
<h2>Calculations Overview</h2>
<p>
    Structural biology analyses begin with the placement of amino acid variants on the individual 3D protein structures
    selected by the pipeline.  From these placements, the pipeline estimates the energetic impact of each amino acid
    substitution via a ΔΔGfolding calculation (<a
        href="https://doi.org/10.1021/acs.jctc.6b00819" target="_blank">Park et al. 2016</a>, <a
        href="https://doi.org/10.3389/fbioe.2020.558247" target="_blank">Frenz et al. 2020</a>).
    Our PathProx algorithm(<a
    href="https://doi.org/10.1016/j.ajhg.2018.01.017" target="_blank">Sivley et al. 2018</a>, <a
        href="https://doi.org/10.1186/s12859-018-2010-z" target="_blank">Sivley et al. 2018</a>) predicts pathogenic variants when
    they better fit with patterns of “known pathogenic” sites vs. random or benign variant sites.   The potential of
    a variant to disrupt protein-protein interaction surfaces(<a
        href="https://doi.org/10.1038/s41592-022-01490-7" target="_blank">Tubiana et al. 2022</a>)
    or perturb post-translational modification sites (<a
        href="https://doi.org/10.1093/bioinformatics/btx496" target="_blank">Wang et al. 2017</a>) is
    also predicted in context of protein 3D structure.
</p><p>
    Genomics is integrated in the pipeline.  Many of our calculations interrogate the Human Genome(<a
        href="https://doi.org/10.1093/nar/gkab1049" target="_blank">Cunningham et al. 2022</a> ) and genomic
    databases.  Mechanically, to map genomic changes to 3D structures, our code navigates technical challenges
    in variant effect prediction and transcript curation(<a
        href="https://doi.org/10.1186/s13059-016-0974-4" target="_blank">McLaren et al. 2016</a>,
        <a href="https://doi.org/10.1093/nar/gkac1052" target="_blank">Bateman et al. 2023</a>).
    Several of our predictive calculations
    integrate sequence constraint, gleaned from both multi-species sequence alignments(<a
    href="https://doi.org/10.1093/bioinformatics/18.suppl_1.S71">Pupko et al. 2002</a>, <a
        href="https://doi.org/10.1093/nar/gkw408" target="_blank">Ashkenazy et al. 2016</a>) and human population
    sequences(<a
    href="https://doi.org/10.1038/s41467-022-30936-x" target="_blank">Li et al. 2022</a>
    ).  Clinvar(<a
            href="https://doi.org/10.1093/nar/gkx1153" target="_isblank">Landrum et al. 2018</a>),
        COSMIC(<a
            href="https://doi.org/10.1093/nar/gky1015" target="_blank>">Tate et al. 2019</a>)
            and gnomAD(<a
            href="https://doi.org/10.1038/s41586-020-2308-7" target="_blank">Karczewski et al. 2020</a>)
    databases are mined for variants needed by
    PathProx’s mathematical spatial analysis, as well as for web-based visualizations.
</p>
    <p>
    Prediction of digenic disease pairs has been an ongoing area of research in both the Capra and Meiler labs.
        In addition to our DiGePred algorithm(17), we also run a newer pair prediction tool (18) against
        each case’s gene list.
    </p>
    <p>
        Please explore the Calculation Details further to learn about the various analyses and their interpretation.
    </p>
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4" crossorigin="anonymous"></script>
<script src="js/includeHTML.js">
</script>

<script>
    includeHTML();
</script>
</div>
</body>
</html>