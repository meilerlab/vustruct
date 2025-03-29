<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>PPI Site Prediction</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">
    <link href="./css/vustruct.css" rel="stylesheet">
</head>
<body>
<div class="container">
    <?php echo file_get_contents("./navbar.html") ?>
<h2>Protein-Protein Interaction (PPI) Site Prediction</h2>
<p>
Many proteins function through binding to other proteins partners, and approximately 60% of disease-associated missense
mutations have been noted to perturb PPIs(<a
        href="https://doi.org/10.1016/j.cell.2015.04.013" target="_blank">Sahni et al. 2015</a>).
Amino acid variants in interaction surfaces could have negligible
impact in folding energetics and still be deleterious through disruption of protein binding to usual partners.
The ScanNet(<a
        href="https://doi.org/10.1038/s41592-022-01490-7" target="_blank">Tubiana et al. 2022</a>) machine learning
algorithm was trained through analysis of the Dockground(<a
        href="https://doi.org/10.1007/978-1-0716-0708-4_17" target="_blank">Kundrotas et al. 2020</a>) database of
3D protein-protein interacting structures.  Operationally, for each variant position covered by an Alphafold model,
we ask ScanNet to predict whether the position participates in a PPI binding surface.
</p><p>
    We report any variant position
to which ScanNet assigns at least 50% probability of being involved in a PPI.
</p>

<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4" crossorigin="anonymous"></script>
<script src="js/includeHTML.js">
</script>

<script>
    includeHTML();
</script>
</body>
</html>