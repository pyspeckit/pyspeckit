<?php
echo "<center>";
$links = array('Home'=>'../index.html', 'Examples'=>'examples.html', 'Documentation'=>'sphinx/index.html');
$leftbrace = "<span class=leftbraces>[</span>";
$rightbrace = "<span class=rightbraces>]</span>";
echo "<nav>";
echo "<ul>";
foreach ($links as $key => $value) {
    echo "<li><a class=navbar href=\"$value\">$key</a></li>";
}
echo "</ul>";
echo "</nav>";

echo "</center>";

echo "<script type=\"text/javascript\">";
echo "<script type=\"text/javascript\">";
echo "";
echo "  var _gaq = _gaq || [];";
echo "  _gaq.push(['_setAccount', 'UA-6253200-7']);";
echo "  _gaq.push(['_trackPageview']);";
echo "";
echo "  (function() {";
echo "    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;";
echo "    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';";
echo "    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);";
echo "  })();";
echo "";
echo "</script>";
?>

