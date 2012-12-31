 // THE LINE PLOT

/*
var data = d3.json( "{{ url_for('serve_data') }}",
	function(i) {return {x: i.x, y: i.y}; });
*/

var data = d3.range(40).map(function(i) {
  return {x: i / 39, y: (Math.sin(i / 3) + 2) / 4};
});


var margin = {top: 10, right: 10, bottom: 20, left: 40},
    width = 600 - margin.left - margin.right,
    height = 300 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .domain([0, 1])
    .range([0, width]);

var y = d3.scale.linear()
    .domain([0, 1])
    .range([height, 0]);

var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y)
    .orient("left");

var line = d3.svg.line()
    .x(function(d) { return x(d.x); })
    .y(function(d) { return y(d.y); });

var svg = d3.select("#main_hero").append("svg")
    .datum(data)
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(0," + height + ")")
    .call(xAxis);

svg.append("g")
    .attr("class", "y axis")
    .call(yAxis);

svg.append("path")
    .attr("class", "line")
    .attr("d", line);