var React = require('react/addons');
var ReactTransitionGroup = React.addons.TransitionGroup;


// CSS
require('normalize.css');
require('../styles/main.css');

require('d3')
var topojson = require('topojson')
var asyncCartogram = require('../../async-cartogram/cartogramaster.js')
var topojsonData = require('json!../../data/arrondissements.json')
var _ = require('underscore.deferred')

var parisData = require('dsv!../../data/metriques-paris.csv')
// Transform this raw csv for the cartogram input
var metrics = Object.keys(parisData[0])
metrics.shift()
var data = {}
metrics.map(function(metric){
  data[metric] = {}
  parisData.map(function(line){
    data[metric][line['Code géographique']] = +line[metric].replace(/,/g, '.')
  })
})
// end transform

//Some vars for our d3 svg map
var scaling = 200000, center = [2.2940220935081865, 48.874100811293694]

var AsyncCartogramApp = React.createClass({
  getInitialState: function(){
    return {metric: null}
  },
  render: function() {
    return (
      <div className='main'>
        <ul id="menu">{metrics.map(m => {
          return <li
                    onMouseOver={this.mouseOverLi.bind(this, m)}
                    onMouseOut={this.mouseOverLi.bind(this, null)}>
                  {m}
                 </li>
        })}</ul>
        <div id="playground">A map will replace me --------> :-o </div>
      </div>
    );
  },

  mouseOverLi: function(m){
    this.setState({metric: m})
  },

  componentDidUpdate: function(){
    this.displayMap()
  },

  componentDidMount: function(){ var _this = this;

    // compute the cartogram geometries for each metric, then cache them.

    var promiseOfGeos = asyncCartogram({
        topology: topojsonData,
        geometries: topojsonData.objects.arrondissements.geometries,
        projection: {
          name: 'mercator', //TODO
          scaling: scaling,
          center: center
        }
      },
      data,
      'c_arinsee'
    );

    // Keep track of the progress (incremented for each task success)
    promiseOfGeos.progress(function(value){
      document.querySelector('#playground').innerHTML = Math.round(value * 100) + "%"
    })

    // All is done
    promiseOfGeos.then(function(a){
      //cache the results
      _this.computedCartograms = a
      _this.displayMap()
    }).fail(function( err ){
      console.log(err.message) // "Oops!"
    })
  },

  displayMap: function(){
    var playground = document.querySelector('#playground')
    playground.innerHTML = ''
    var svg = d3.select(playground).append("svg").attr("id", "leSVG")
    var path, featureCollection;
    if (this.state.metric == null){// Show the real map

      var projection = d3.geo.mercator()
      .center(center)
      .scale(scaling)

      path = d3.geo.path().projection(projection)
      //convert to geojson
      featureCollection = topojson.feature(topojsonData, topojsonData.objects.arrondissements);
    } else { //Draw the population cartogram
      var a = this.computedCartograms
      featureCollection = a[this.state.metric] //get first task result

      // path with identity projection
      path = d3.geo.path()
      .projection(null);
    }

    // Now display the cartogram

    var nodes = []

    //computing the area centroids from their geometry
    featureCollection.features.forEach((d, i) => {
      var centroid = path.centroid(d);
      if (centroid.some(isNaN)) {
        return;
      }
      centroid.x = centroid[0];
      centroid.y = centroid[1];
      centroid.feature = d;
      nodes.push(centroid);
    });

    //adding the map areas
    var areas = svg.append('g').attr('class', 'nodes').selectAll("g")
      .data(nodes)
      .enter().append("g")
    areas
      .append("path")
      .attr("d", function(d) { return path(d.feature); })


    //let's add names to our map areas
    areas
      .append("text")
      .attr("x", d => d.x)
      .attr("y", d => d.y)
      .text(d => d.feature.properties.c_ar)

    //let's add la Tour Eiffel
    areas.each(function(d){
      if (d.feature.properties.c_ar == 7){
        d3.select(this)
            .append("image")
            .attr({
              "x":  d => d.x - 35,
              "y":  d => d.y - 35,
              "width":  "20px",
              "height":  "20px",
              "xlink:href": "../images/latoureiffel.svg"
            })
      }
    })
  }
});
React.render(<AsyncCartogramApp />, document.getElementById('content')); // jshint ignore:line

module.exports = AsyncCartogramApp;
