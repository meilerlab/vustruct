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

<h2>Patient Case Tutorial - VCF input</h2>
    <p>Variant Call Format (VCF) is the standard tab-delimited file format for communicating alternate alleles which
    differ from the reference genome</p>
    <p>In the simple example below.......  </p>

    <pre><code>
$ cat RTEL1.vcf
#CHROM	POS	ID	REF	ALT
20	63659451	.	C	T	.	.	.
20	63662544	.	A	G	.	.	.
20	63667545	.	G	T	.	.	.
20	63687668	.	C	T	.	.	.
20	63687765	.	G	T	.	.	.
20	63688001	.	G	C	.	.	.
20	63688578	.	G	T	.	.	.
20	63689115	.	G	A	.	.	.
20	63689583	.	C	A	.	.	.
20	63689821	.	C	G	.	.	.
20	63694906	.	T	C	.	.	.
20	63695093	.	A	C	.	.	.
20	63695619	.	G	A	.	.	.

    </pre></code>
    <p>

        To Do:  Down a file link, instructions, screen-shots.

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