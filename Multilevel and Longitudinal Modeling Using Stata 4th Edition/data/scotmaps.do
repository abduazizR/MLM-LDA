** corrections from German Rodriguez

// Do-file to create the choropleth maps shown in figures 6.1 and 6.2 of
// Multilevel and Longitudinal Modeling Using Stata, Sophia Rabe-Hesketh and
// Anders Skrondal, Stata Press, 2005.

// Use the dataset containing the polygon outlines of Scottish counties.

use https://www.stata-press.com/data/mlmus/scotmap

// Merge the county outlines to a dataset of various SMR estimates by county.

merge regid using https://www.stata-press.com/data/mlmus/smrs


// For two of the SMR estimates, crude and normal, code the SMRs into
// variables that measure the SMR in increasing levels from 1 to 5.

gen  crude_level = irecode(smr_crude,  -1, 50, 80, 120, 200, 700)
gen normal_level = irecode(smr_normal, -1, 50, 80, 120, 200, 700)

// Create a choropleth map where each county's fill color (in this
// case level of gray) is determined by its crude SMR level.  

// We draw 5 plots of county outlines, one for each crude SMR level,
// and fill the outlines of each plot with a different shade of gray,
// e.g., option fc(gs12) is an abbreviation of fcolor(gs12) and selects
// a light gray fill color.  The lc(gs0) options are short for
// lcolor(gs0) and select black outlines for all regions and lw(vthin)
// is short for lwidth(vthin) and selects a very thin outline for the
// regions.  The option nodrop is short for -nodropbase- meaning
// connect the regions as stand-alone polygons, rather than a level
// above the x-asis of the graph.

// Other options ensure the proper aspect ratio for the map, turn off default
// axes, set the graph size, and position and annotate the legend.  See 
// -help twoway options- for details.

local options nodropb lc(gs0) lw(vthin) cmissing(n)
twoway                                                     ///
   area lat lon if crude_level == 1, fc(gs12) `options' || ///
   area lat lon if crude_level == 2, fc(gs9)  `options' || ///
   area lat lon if crude_level == 3, fc(gs6)  `options' || ///
   area lat lon if crude_level == 4, fc(gs3)  `options' || ///
   area lat lon if crude_level == 5, fc(gs0)  `options'    ///
      aspect(1.6) scheme(s1mono) xscale(off) yscale(off)           ///
      plotregion(lcolor(none) margin(l=4)) ysize(7) xsize(5.5)     ///
      legend(position(11) ring(0) cols(1) size(small) symxsize(4)  ///
      order(5 "200 to 700" 4 "120 to 200" 3 "80 to 120"            ///
            2 "50 to 80"   1 "0 to 50"))

*twoway                                                                      ///
*   area lat lon if crude_level == 1, nodropb fc(gs12) lc(gs0) lw(vthin) ||  ///
*   area lat lon if crude_level == 2, nodropb fc(gs9)  lc(gs0) lw(vthin) ||  ///
*   area lat lon if crude_level == 3, nodropb fc(gs6)  lc(gs0) lw(vthin) ||  ///
*   area lat lon if crude_level == 4, nodropb fc(gs3)  lc(gs0) lw(vthin) ||  ///
*   area lat lon if crude_level == 5, nodropb fc(gs0)  lc(gs0) lw(vthin)     ///
*        aspect(1.6) scheme(s1mono) xscale(off) yscale(off)                  ///
*        plotregion(lcolor(none) margin(l=4)) ysize(7) xsize(5.5)            ///
*        legend(position(11) ring(0) cols(1) size(small) symxsize(4)         ///
*               order(5 "200 to 700" 4 "120 to 200" 3 "80 to 120"            ///
*                     2 "50 to 80"   1 "0 to 50"))

// Create a similar choropleth map where the fill colors are determined by
// each county's normal estimate of SMR level.

local options nodropb lc(gs0) lw(vthin) cmissing(n)
twoway                                                     ///
   area lat lon if normal_level == 1, fc(gs12) `options' || ///
   area lat lon if normal_level == 2, fc(gs9)  `options' || ///
   area lat lon if normal_level == 3, fc(gs6)  `options' || ///
   area lat lon if normal_level == 4, fc(gs3)  `options' || ///
   area lat lon if normal_level == 5, fc(gs0)  `options'    ///
      aspect(1.6) scheme(s1mono) xscale(off) yscale(off)           ///
      plotregion(lcolor(none) margin(l=4)) ysize(7) xsize(5.5)     ///
      legend(position(11) ring(0) cols(1) size(small) symxsize(4)  ///
      order(5 "200 to 700" 4 "120 to 200" 3 "80 to 120"            ///
            2 "50 to 80"   1 "0 to 50"))
            
*twoway                                                                      ///
*   area lat lon if normal_level == 1, nodropb fc(gs12) lc(gs0) lw(vthin) || ///
*   area lat lon if normal_level == 2, nodropb fc(gs9)  lc(gs0) lw(vthin) || ///
*   area lat lon if normal_level == 3, nodropb fc(gs6)  lc(gs0) lw(vthin) || ///
*   area lat lon if normal_level == 4, nodropb fc(gs3)  lc(gs0) lw(vthin) || ///
*   area lat lon if normal_level == 5, nodropb fc(gs0)  lc(gs0) lw(vthin)    ///
*        aspect(1.6) scheme(s1mono) xscale(off) yscale(off)                  ///
*        plotregion(lcolor(none) margin(l=4)) ysize(7) xsize(5.5)            ///
*        legend(position(11) ring(0) cols(1) size(small) symxsize(4)         ///
*               order(5 "200 to 700" 4 "120 to 200" 3 "80 to 120"            ///
*                     2 "50 to 80"   1 "0 to 50"))

