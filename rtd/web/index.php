<?php // Page Template version 3.3

// PLEASE UPDATE THESE VALUES FOR EVERY PAGE
// -----------------------------------------------------------------------------
$title        = 'GHOST Data Reduction Software Documentation Index';
$description  = 'The index page for the GHOST Data Reduction Software Documentation';
$subject      = 'GHOST';

// Include configuration file & begin the XHTML response
// -----------------------------------------------------------------------------
include ($_SERVER['DOCUMENT_ROOT'] .'/rsaa/config.php');

echo $DocType;
echo $Head;
// insert additional header statements here //
echo $Body;
echo $Banner;
//include $Menu;//
include ($_SERVER['DOCUMENT_ROOT'] .$Menu);

echo '<style>#body {margin-left:0}</style>'


// BEGIN: Page Content
// =============================================================================
?>
<!-- START MAIN PAGE CONTENT -->
<img src="GHOST_logos.jpg">
<?php
include 'index.html'
?>
<!-- END MAIN PAGE CONTENT -->
<?php
// =============================================================================
// END: Page Content


// Complete the XHTML response
// -----------------------------------------------------------------------------
echo $Update;
echo $Analytics;
echo $Footer;
?>

