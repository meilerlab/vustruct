<!doctype html>
<html lang="en">

<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>VUStruct Documentation</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">
    <link href="./css/vustruct.css" rel="stylesheet">
</head>
<body>
<div class="container">

<?php echo file_get_contents("./navbar.html") ?>

<h1>Welcome to VUStruct</h1>
    <p>

        Mutations in protein coding regions can be deleterious through various
        mechanisms(<a href="https://doi.org/10.1139/bcb-2018-0007" target="_blank">Taipale 2019<a>).  To holistically
        inform structural biologists, VUStruct accepts sets of Variants of Unknown Significance (VUSs) as either
        a VCF file of human genomic coordinates or as a CSV list of amino acid changes to specific Uniprot proteins.
        The pipeline automatically selects protein structures which cover missense variants, launches the compute
        jobs on Vanderbilt's <a href="https://www.vanderbilt.edu/accre/" target="_blank">ACCRE</a> computer cluster,
        and creates a final drill-downable <a href="https://staging.meilerlab.org/vustruct/SampleCase/">report page</a>

    </p>
    <p>
    These documentation pages provide pipeline input examples, output report examples, and documents our thinking
    on how the suite of computational methods, together with the generated website, inform hypothesis generation.
</p>
    <p>
    <strong>
    The VUStruct compute pipeline is intended only for academic research, and is provided free of charge to all users, 
    including commercial users.</strong></p>
</p>

<div class="row text-center">

<image src="./VUStructSchematic.png" class="img-fluid align-content-center"  alt="VUStruct Pipeline Schematic" style="padding: 1em"></image>
</div>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4" crossorigin="anonymous"></script>
<script src="js/includeHTML.js">
</script>

<script>
    includeHTML();
</script>

</div>
</body>
</html>
