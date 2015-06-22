var React = require('react/addons');
var ReactTransitionGroup = React.addons.TransitionGroup;

require('d3')
var topojson = require('topojson')
var asyncCartogram = require('../../async-cartogram/cartogramaster.js')
var topojsonData = require('json!../../data/arrondissements.json')
var _ = require('underscore.deferred')
var parisData = require('dsv!../../data/metriques-paris.csv')

var metrics = Object.keys(parisData[0])
metrics.shift()
var data = {}
metrics.map(function(metric){
  data[metric] = {}
  parisData.map(function(line){
    data[metric][line['Code géographique']] = +line[metric].replace(/,/g, '.')
  })
})
debugger;

// CSS
require('normalize.css');
require('../styles/main.css');

var AsyncCartogramApp = React.createClass({
  getInitialState: function(){
    return {
      trueMap: false
    }
  },
  render: function() {
    return (
      <div className='main'>
        <ul id="menu">{metrics.map(m => <li>{m}</li>)}</ul>
        <div id="playground">A map will replace me --------> :-o </div>
      </div>
    );
  },

  componentDidMount: function(){ var _this = this;

    var scaling = 200000, center = [2.2940220935081865, 48.874100811293694]

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

    promiseOfGeos.progress(function(value){
      document.querySelector('#playground').innerHTML = Math.round(value * 100) + "%"
    })

    promiseOfGeos.then(function(a){
      var playground = document.querySelector('#playground')
      playground.innerHTML = ''
      var svg = d3.select(playground).append("svg").attr("id", "leSVG")

      var path, featureCollection;
      debugger;
      if (_this.state.trueMap){// Show the real geographical map

        var projection = d3.geo.mercator()
        .center(center)
        .scale(scaling)

        //var centerPoint = projection([2.2940220935081865, 48.874100811293694])

        path = d3.geo.path().projection(projection)
        //convert to geojson
        featureCollection = topojson.feature(topojsonData, topojsonData.objects.arrondissements);
      } else { //Draw the population cartogram

        featureCollection = a[Object.keys(a)[0]] //get first task result

        // path with identity projection
        path = d3.geo.path()
        .projection(null);
      }

      // Now display the cartogram

      var nodes = []

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

      var nodeGroup = svg.append('g').attr('class', 'nodes')

      var node = nodeGroup.selectAll("g")
      .data(nodes)
      .enter().append("g")
      .attr("transform", function(d) { return "translate(" + -d.x + "," + -d.y + ")"; })
      .append("path")
      .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
      .attr("d", function(d) { return path(d.feature); })


    }).fail(function( err ){
      console.log(err.message); // "Oops!"
    });
  }
});
React.render(<AsyncCartogramApp />, document.getElementById('content')); // jshint ignore:line

module.exports = AsyncCartogramApp;
