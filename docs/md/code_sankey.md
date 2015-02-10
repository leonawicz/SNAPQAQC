# Code Flow



## Related items
SNAP QAQC code feeds into many other projects, including the SNAP data QAQC master R Shiny app, Community Charts apps, and ALFRESCO-related projects.
The primary order QAQC code is most directly related to data extraction.
Secondary order code relates to organization and preparation of useful data products which are taken up by web applications and other projects to investigate the data.
Tertiary order code relates to specific QAQC analyses.

### Files and Data
This project makes use of any SNAP data, particularly gridded climate, fire, and vegetation data sets.
Files generated as part of SNAP QAQC include fortified datasets which are more accessible, compact, easier to manipulate, understand, and analyze for specific purposes.

### SNAP QAQC code flow











<iframe srcdoc=' &lt;!doctype HTML&gt;
&lt;meta charset = &#039;utf-8&#039;&gt;
&lt;html&gt;
  &lt;head&gt;
    &lt;link rel=&#039;stylesheet&#039; href=&#039;http://timelyportfolio.github.io/rCharts_d3_sankey/css/sankey.css&#039;&gt;
    
    &lt;script src=&#039;http://d3js.org/d3.v3.min.js&#039; type=&#039;text/javascript&#039;&gt;&lt;/script&gt;
    &lt;script src=&#039;http://timelyportfolio.github.io/rCharts_d3_sankey/js/sankey.js&#039; type=&#039;text/javascript&#039;&gt;&lt;/script&gt;
    
    &lt;style&gt;
    .rChart {
      display: block;
      margin-left: auto; 
      margin-right: auto;
      width: 900px;
      height: 800px;
    }  
    &lt;/style&gt;
    
  &lt;/head&gt;
  &lt;body &gt;
    
    &lt;div id = &#039;chart818414c563b&#039; class = &#039;rChart d3_sankey&#039;&gt;&lt;/div&gt;    
    ï»¿&lt;!--Attribution:
Mike Bostock https://github.com/d3/d3-plugins/tree/master/sankey
Mike Bostock http://bost.ocks.org/mike/sankey/
--&gt;

&lt;script&gt;
(function(){
var params = {
 &quot;dom&quot;: &quot;chart818414c563b&quot;,
&quot;width&quot;:    900,
&quot;height&quot;:    800,
&quot;data&quot;: {
 &quot;source&quot;: [ &quot;AR4_AR5_extract.slurm&quot;, &quot;CRU_extract.slurm&quot;, &quot;AR4_AR5_extract.R&quot;, &quot;AR4_AR5_extract.R&quot;, &quot;AR4_AR5_extract.R&quot;, &quot;CRU_extract.R&quot;, &quot;CRU_extract.R&quot;, &quot;CRU_extract.R&quot;, &quot;shapes2cells_AKCAN2km_5pct.RData&quot;, &quot;Community location data: locs.csv&quot;, &quot;2-km AK-CAN downscaled climate data&quot;, &quot;shapes2cells_AKCAN2km_5pct.RData&quot;, &quot;Community location data: locs.csv&quot;, &quot;2-km AK-CAN downscaled climate data&quot;, &quot;Regional shapefiles&quot;, &quot;shapes2cells.R&quot;, &quot;shapes2cells.R&quot;, &quot;shapes2cells.R&quot;, &quot;cities_setup.slurm&quot;, &quot;stats_setup.R&quot;, &quot;samples_setup.R&quot;, &quot;samples_setup_CRU31.R&quot;, &quot;qaqc_app_metadata.R&quot;, &quot;qaqc_app_metadata.R&quot;, &quot;cities_meta_akcan2km.RData&quot;, &quot;alfStatsByRep_Rmpi.slurm&quot;, &quot;shapes2cells_AKCAN1km.RData&quot;, &quot;shapes2cells_AKCAN1km_rmNA.RData&quot;, &quot;alfStatsByRep_Rmpi.R&quot;, &quot;alfStatsByRep.R&quot;, &quot;Intermediary data&quot;, &quot;getAlfStatsAndDensities.R&quot;, &quot;getAlfStatsAndDensities.R&quot;, &quot;getAlfStatsAndDensities.R&quot;, &quot;1-km ALFRESCO outputs&quot;, &quot;1-km ALFRESCO outputs&quot; ],
&quot;target&quot;: [ &quot;AR4_AR5_extract.R&quot;, &quot;CRU_extract.R&quot;, &quot;cities_setup.R&quot;, &quot;stats_setup.R&quot;, &quot;samples_setup.R&quot;, &quot;cities_setup_CRU31.R&quot;, &quot;stats_setup_CRU31.R&quot;, &quot;samples_setup_CRU31.R&quot;, &quot;AR4_AR5_extract.R&quot;, &quot;AR4_AR5_extract.R&quot;, &quot;AR4_AR5_extract.R&quot;, &quot;CRU_extract.R&quot;, &quot;CRU_extract.R&quot;, &quot;CRU_extract.R&quot;, &quot;shapes2cells.R&quot;, &quot;shapes2cells_AKCAN2km_5pct.RData&quot;, &quot;shapes2cells_AKCAN1km.RData&quot;, &quot;shapes2cells_AKCAN1km_rmNA.RData&quot;, &quot;cities_setup.R&quot;, &quot;Master QAQC app: meta.RData&quot;, &quot;Master QAQC app: meta.RData&quot;, &quot;Master QAQC app: meta.RData&quot;, &quot;Master QAQC app: meta.RData&quot;, &quot;cities_meta_akcan2km.RData&quot;, &quot;Master QAQC app: meta.RData&quot;, &quot;alfStatsByRep_Rmpi.R&quot;, &quot;alfStatsByRep_Rmpi.R&quot;, &quot;alfStatsByRep_Rmpi.R&quot;, &quot;alfStatsByRep.R&quot;, &quot;Intermediary data&quot;, &quot;getAlfStatsAndDensities.R&quot;, &quot;Master QAQC app: meta.RData&quot;, &quot;App data&quot;, &quot;Other QAQC analyses&quot;, &quot;alfStatsByRep_Rmpi.R&quot;, &quot;alfStatsByRep.R&quot; ],
&quot;value&quot;: [      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1 ] 
},
&quot;nodeWidth&quot;:     15,
&quot;nodePadding&quot;:     10,
&quot;layout&quot;:     32,
&quot;margin&quot;: {
 &quot;right&quot;:     20,
&quot;left&quot;:     20,
&quot;bottom&quot;:     50,
&quot;top&quot;:     50 
},
&quot;title&quot;: &quot;Code Flow&quot;,
&quot;id&quot;: &quot;chart818414c563b&quot; 
};

params.units ? units = &quot; &quot; + params.units : units = &quot;&quot;;

//hard code these now but eventually make available
var formatNumber = d3.format(&quot;0,.0f&quot;),    // zero decimal places
    format = function(d) { return formatNumber(d) + units; },
    color = d3.scale.category20();

if(params.labelFormat){
  formatNumber = d3.format(&quot;.2%&quot;);
}

var svg = d3.select(&#039;#&#039; + params.id).append(&quot;svg&quot;)
    .attr(&quot;width&quot;, params.width)
    .attr(&quot;height&quot;, params.height);
    
var sankey = d3.sankey()
    .nodeWidth(params.nodeWidth)
    .nodePadding(params.nodePadding)
    .layout(params.layout)
    .size([params.width,params.height]);
    
var path = sankey.link();
    
var data = params.data,
    links = [],
    nodes = [];
    
//get all source and target into nodes
//will reduce to unique in the next step
//also get links in object form
data.source.forEach(function (d, i) {
    nodes.push({ &quot;name&quot;: data.source[i] });
    nodes.push({ &quot;name&quot;: data.target[i] });
    links.push({ &quot;source&quot;: data.source[i], &quot;target&quot;: data.target[i], &quot;value&quot;: +data.value[i] });
}); 

//now get nodes based on links data
//thanks Mike Bostock https://groups.google.com/d/msg/d3-js/pl297cFtIQk/Eso4q_eBu1IJ
//this handy little function returns only the distinct / unique nodes
nodes = d3.keys(d3.nest()
                .key(function (d) { return d.name; })
                .map(nodes));

//it appears d3 with force layout wants a numeric source and target
//so loop through each link replacing the text with its index from node
links.forEach(function (d, i) {
    links[i].source = nodes.indexOf(links[i].source);
    links[i].target = nodes.indexOf(links[i].target);
});

//now loop through each nodes to make nodes an array of objects rather than an array of strings
nodes.forEach(function (d, i) {
    nodes[i] = { &quot;name&quot;: d };
});

sankey
  .nodes(nodes)
  .links(links)
  .layout(params.layout);
  
var link = svg.append(&quot;g&quot;).selectAll(&quot;.link&quot;)
  .data(links)
.enter().append(&quot;path&quot;)
  .attr(&quot;class&quot;, &quot;link&quot;)
  .attr(&quot;d&quot;, path)
  .style(&quot;stroke-width&quot;, function (d) { return Math.max(1, d.dy); })
  .sort(function (a, b) { return b.dy - a.dy; });

link.append(&quot;title&quot;)
  .text(function (d) { return d.source.name + &quot; â†’ &quot; + d.target.name + &quot;\n&quot; + format(d.value); });

var node = svg.append(&quot;g&quot;).selectAll(&quot;.node&quot;)
  .data(nodes)
.enter().append(&quot;g&quot;)
  .attr(&quot;class&quot;, &quot;node&quot;)
  .attr(&quot;transform&quot;, function (d) { return &quot;translate(&quot; + d.x + &quot;,&quot; + d.y + &quot;)&quot;; })
.call(d3.behavior.drag()
  .origin(function (d) { return d; })
  .on(&quot;dragstart&quot;, function () { this.parentNode.appendChild(this); })
  .on(&quot;drag&quot;, dragmove));

node.append(&quot;rect&quot;)
  .attr(&quot;height&quot;, function (d) { return d.dy; })
  .attr(&quot;width&quot;, sankey.nodeWidth())
  .style(&quot;fill&quot;, function (d) { return d.color = color(d.name.replace(/ .*/, &quot;&quot;)); })
  .style(&quot;stroke&quot;, function (d) { return d3.rgb(d.color).darker(2); })
.append(&quot;title&quot;)
  .text(function (d) { return d.name + &quot;\n&quot; + format(d.value); });

node.append(&quot;text&quot;)
  .attr(&quot;x&quot;, -6)
  .attr(&quot;y&quot;, function (d) { return d.dy / 2; })
  .attr(&quot;dy&quot;, &quot;.35em&quot;)
  .attr(&quot;text-anchor&quot;, &quot;end&quot;)
  .attr(&quot;transform&quot;, null)
  .text(function (d) { return d.name; })
.filter(function (d) { return d.x &lt; params.width / 2; })
  .attr(&quot;x&quot;, 6 + sankey.nodeWidth())
  .attr(&quot;text-anchor&quot;, &quot;start&quot;);

// the function for moving the nodes
  function dragmove(d) {
    d3.select(this).attr(&quot;transform&quot;, 
        &quot;translate(&quot; + (
                   d.x = Math.max(0, Math.min(params.width - d.dx, d3.event.x))
                ) + &quot;,&quot; + (
                   d.y = Math.max(0, Math.min(params.height - d.dy, d3.event.y))
                ) + &quot;)&quot;);
        sankey.relayout();
        link.attr(&quot;d&quot;, path);
  }
})();
&lt;/script&gt;
    
    
    &lt;script&gt;
      var cscale = d3.scale.category20b();
      d3.selectAll(&#039;#chart818414c563b svg path.link&#039;)
        .style(&#039;stroke&#039;, function(d){
          return cscale(d.source.name);
        })
      d3.selectAll(&#039;#chart818414c563b svg .node rect&#039;)
        .style(&#039;fill&#039;, function(d){
          return cscale(d.name)
        })
        .style(&#039;stroke&#039;, &#039;none&#039;)
    &lt;/script&gt;
        
  &lt;/body&gt;
&lt;/html&gt; ' scrolling='no' frameBorder='0' seamless class='rChart  http://timelyportfolio.github.io/rCharts_d3_sankey/libraries/widgets/d3_sankey  ' id='iframe-chart818414c563b'> </iframe>
 <style>iframe.rChart{ width: 100%; height: 400px;}</style>
<style>iframe.rChart{ width: 100%; height: 840px;}</style>
