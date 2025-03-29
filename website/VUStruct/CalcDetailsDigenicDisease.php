<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Digenic Disease Prediction</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">
    <link href="./css/vustruct.css" rel="stylesheet">
</head>
<body>
<div class="container">
    <?php echo file_get_contents("./navbar.html") ?>
<h2>DiGePred and DIEP</h2>
<p>
Disease emergence through compound heterozygous variation is well understood.
Individually, a Proband’s disease-free parents’ gene functions enough under single variation.
But the combination of gene variants from each parent triggers creates the observed phenotype in the proband.
Thus, compound heterozygous variant pairs are always given special attention in each week’s UDN review.
</p>
<p>
Multigenic disease can analogously arise when two genes in a single pathway are damaged via variations inherited from each parent.
</p>
    <p>
An increasing database of annotated digenic(<a
    href="https://doi.org/10.1093/nar/gkv1068" target="_blank">Gazzo et al. 2016</a>) and multigenic(<a
    href="https://doi.org/10.1093/database/baac023" target="_blank">Nachtegael et al. 2022</a>)
diseases is being curated from the literature.
Former Ph.D. student Souhrid Mukherjee asked whether latent digenic diseases might be predicted through machine
learning strategies, and he produced a tool called DiGePred(<a
        href="https://doi.org/10.1016/j.ajhg.2021.08.010" target="_blank">Mukherjee et al. 2021</a>) Another group
    broadened the training inputs, while retaining a similar machine architecture, to create DIEP(<a
        href="https://doi.org/10.1016/j.csbj.2022.07.011" target="_blank">Yuan et al. 2022</a>)
    </p>

    <p>Add sample figure here.</p>

<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4" crossorigin="anonymous"></script>
<script src="js/includeHTML.js">
</script>

<script>
    includeHTML();
</script>
</body>
</html>
