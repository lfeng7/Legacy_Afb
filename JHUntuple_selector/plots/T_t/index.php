<//-- 
author currently uknown although my suspicion is 
that this was written by Giovanni Petrucciani

Code received from N. Tran, who received it from 
Cristina Botta

   To whoever wrote this: thanks!!!!!!
--!>

<html>
<head>
<title><?php echo getcwd(); ?></title>
<style type='text/css'>
   div.pic { 
   display: block;
   float: left;
   border: 0px solid gray;
   margin: 3px;
 }
</style>
</head>
<body>
<!-- HERE COMMENTED REGION BEGINS... 
<h1><?php echo getcwd(); ?></h1>
... HERE COMMENTED REGION ENDS -->
<h2><a name="plots">Plots</a></h2>

<p><form>Filter: <input type="text" name="match" size="30" value="<?php if (isset($_GET['match'])) print htmlspecialchars($_GET['match']);  ?>" /><input type="Submit" value="Go" /></form></p>

<div>
<?php
$remainder = array();
if ($_GET['noplots']) {
  print "Plots will not be displayed.\n";
} else {
  //$other_exts = array('.pdf', '.cxx', '.eps');

  print $myTest;

  foreach (glob("*") as $filename) {
    if (!preg_match('/.*\.png$/', $filename) && !preg_match('/.*\.gif.*$/', $filename)) {
      if (!preg_match('/.*\.php.*$/', $filename)) {
	array_push($remainder, $filename);
      }
      continue;
    }
    if (isset($_GET['match']) && !fnmatch('*'.$_GET['match'].'*', $filename)) continue;

    $basefilename=substr($filename,0,-4);
     print "<div class='pic'>\n";
    print "<h6 align=center> $basefilename <a href=$filename>[.png]</a>";
    
    $testfile=$basefilename.".eps";
    if (file_exists($testfile)){
      print "<a href=\"$basefilename.eps\"> [.eps] </a>";
    }
    $testfile=$basefilename.".root";
    if (file_exists($testfile)){
      print "<a href=\"$basefilename.root\"> [.root] </a>";
    }
    $testfile=$basefilename.".pdf";
    if (file_exists($testfile)){
      print "<a href=\"$basefilename.pdf\"> [.pdf] </a>";
    }
    $testfile=$basefilename.".C";
    if (file_exists($testfile)){
      print "<a href=\"$basefilename.C\"> [.C] </a>";
    }

    print " </h6>";

    print "<img src=\"$filename\" style=\"border: none; width: 300px; \"></a>";
    $others = array();
    foreach ($other_exts as $ex) {
        if (file_exists(str_replace('.png', $ex, $filename))) {
            array_push($others, "<a class=\"file\" href=\"".str_replace('.png', $ex, $filename)."\">[" . $ex . "]</a>");
        }
    }
    if ($others) print "<p>Also as ".implode(', ',$others)."</p>";
    print "</div>";
  }

  print "</table> \n " ;
 }
?>
</div>
<div style="display: block; clear:both;">
  <h2><a name="files">&nbsp;<P>Files</a></h2>
<ul>
<?
foreach ($remainder as $filename) {
  print "<li><a href=\"$filename\">$filename</a></li>";
}
?>
</ul>
</div>
<div>
  <p>
    <?
    foreach(glob("*.txt") as $textfile){
      print "<HR><h6>$textfile</h6><HR><PRE>";
      print file_get_contents($textfile);
      print "</PRE>";
    }
    ?>
  </p>
</div>
</body>
</html>
