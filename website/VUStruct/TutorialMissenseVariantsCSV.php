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

<h2>Missense Variant .csv format Tutorial</h2>
    <p></p>
    <p>In the simple example below.......  </p>

    <pre><code>
$ cat UDN166567_development_missense.csv
,gene,transcript,unp,refseq,mutation
0,AMPD1,ENST00000520113.7,P23109-1,NM_000036.2,K287I
2,CFAP54,ENST00000524981.9,Q96N23-1,NM_001306084.1,F1728S
3,CFAP54,ENST00000524981.9,Q96N23-1,NM_001306084.1,T405A
4,SCN10A,ENST00000449082.3,Q9Y5Y9,NM_001293306.2,D826Y
5,SCN10A,ENST00000449082.3,Q9Y5Y9,NM_001293306.2,P477L
6,COL6A3,ENST00000295550.9,P12111-1,NM_004369.3,G2297A
9,DYNC2H1,ENST00000375735.7,Q8NCM8-1,NM_001377.2,G1222E
11,PRR14L,chr22,31715471,G/A,ENST00000327423.11,Q5THK1-1,NM_173566.2,R790C
12,SAMD9,chr7,93103844,A/T,ENST00000379958.3;ENST00000620985.4,Q5K651,NM_001193307.1,W752R
13,TNFRSF13B,chr17,16940415,G/T,ENST00000261652.7,O14836-1,NM_012452.2,A181E
14,TNFRSF13B,chr17,16940415,G/T,ENST00000583789.1,O14836-2,NA,A135E
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