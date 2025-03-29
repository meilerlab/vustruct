<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>PTM Prediction</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-rbsA2VBKQhggwzxH7pPCaAqO46MgnOM80zW1RWuH61DGLwZJEdK2Kadq2F9CUG65" crossorigin="anonymous">
    <link href="./css/vustruct.css" rel="stylesheet">
</head>
<body>
<div class="container">
    <?php echo file_get_contents("./navbar.html") ?>
<h2>Post-translational Modifications (PTMs)</h2>
<p>
    Post-translational modifications to proteins are often critical to function.
    However, these modifications (phosphorylation, glycosylation, and others) are typically unresolved in
    experimental X-ray structures, and never seen in model structures.  Sequence-based prediction of
    post-translational modification sites has been evolving in recent decades, and we have integrated
    MutsiteDeep (<a
        href="https://doi.org/10.1093/bioinformatics/btx496" target="_blank">Wang et al. 2017</a>), the first deep-learning
    framework for prediction of these sites.
    </p><p>
    We report variants at sites predicted to have at least a 50% likelihood of post-translational modification.
</p>


    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-kenU1KFdBIe4zVF0s0G1M5b4hcpxyD9F7jL+jjXkk+Q2h455rYXK/7HAuoJl+0I4" crossorigin="anonymous"></script>
<script src="js/includeHTML.js">
</script>

<script>
    includeHTML();
</script>
</body>
</html>