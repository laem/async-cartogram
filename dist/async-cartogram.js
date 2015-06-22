(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory();
	else if(typeof define === 'function' && define.amd)
		define(factory);
	else if(typeof exports === 'object')
		exports["AsyncCartogram"] = factory();
	else
		root["AsyncCartogram"] = factory();
})(this, function() {
return /******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};

/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {

/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId])
/******/ 			return installedModules[moduleId].exports;

/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			exports: {},
/******/ 			id: moduleId,
/******/ 			loaded: false
/******/ 		};

/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);

/******/ 		// Flag the module as loaded
/******/ 		module.loaded = true;

/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}


/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;

/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;

/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";

/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(0);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ function(module, exports, __webpack_require__) {

	  /*
	  Cartogram forked from the original http://prag.ma/code/d3-cartogram/

	  - Designed to compute a series of cartograms (called tasks),
	  using web workers, or in a single node thread.

	  - For this, only a subset of d3 is used, and the usage is more constrained.

	  - Also attempts to add the effective range tip from Algorithms for Cartogram Animation,
	  by Ouyang et al., which should make it slightly faster.

	  Usage :
	  --------

	  var promiseOfGeos = cartogramaster(
	    {
	      topology: topojsonData,
	      // the geometries of the GeometryCollection object to reshape :
	      geometries: topojsonData.objects.OBJECTNAME.geometries,
	      projection: {
	        name: 'mercator',
	        scaling: scaling,
	        translation: [X,Y], // optional
	        center: center // optional
	      }
	    },
	    values, // { taskId => { geoJsonFeatureValue => area wanted} }
	    featureProperty // geoJsonFeatureIdKey to link geojson features to values
	  );


	 **** Original cartogram.js usage *********
	 *
	 * var cartogram = d3.cartogram()
	 *  .projection(d3.geo.albersUsa())
	 *  .value(function(d) {
	 *    return Math.random() * 100;
	 *  });
	 * d3.json("path/to/topology.json", function(topology) {
	*  var features = cartogram(topology, topology.objects.OBJECTNAME.geometries);
	 *  d3.select("svg").selectAll("path")
	 *    .data(features)
	 *    .enter()
	 *    .append("path")
	 *      .attr("d", cartogram.path);
	 * });
	 */

	function spawnWorker(){
	  var Wok
	  if (typeof window !== 'undefined'){//browser
	    Wok = __webpack_require__(1);
	  } else {//node
	    Wok = __webpack_require__(3)
	  }
	  return new Wok()
	}


	var _ = __webpack_require__(6)

	function cartogramaster(geo, values, featureProperty) {
	  var dfd = new _.Deferred()

	  // How many cartograms to compute ?
	  var tasks = Object.keys(values)

	  var n = tasks.length

	  var workers = [],
	  results = {};

	  function post(worker, index){
	    worker.postMessage({
	      do: 'carto',
	      geo: geo,
	      values: values[index],
	      featureProperty: featureProperty,
	      task: index
	    })
	  }

	  // Spawn 8 workers max, feed them sequentially until all tasks are done.
	  while (workers.length < 9 && tasks.length > 0){

	   var worker = spawnWorker()
	   workers.push(worker)

	   work(worker, tasks.pop())

	  }

	  function work(worker, i){

	   worker.onmessage = function(event){
	     var data = event
	     if (typeof event.data!== 'undefined') data = event.data

	     if (data.done === 'processing'){
	       dfd.notify(Object.keys(results).length / n)
	       results[data.task] = data

	       //the end
	       if (Object.keys(results).length === n){
	         dfd.resolve(results);
	         workers.forEach(function(w){
	           w.terminate()
	         });
	         return;
	       }
	       if (tasks.length > 0){
	         var j = tasks.pop()
	         post(worker, j)
	       }
	     }
	   }

	   post(worker, i)
	  }

	return dfd
	};


	module.exports = cartogramaster


/***/ },
/* 1 */
/***/ function(module, exports, __webpack_require__) {

	module.exports = function() {
		return __webpack_require__(2)("/******/ (function(modules) { // webpackBootstrap\n/******/ \t// The module cache\n/******/ \tvar installedModules = {};\n\n/******/ \t// The require function\n/******/ \tfunction __webpack_require__(moduleId) {\n\n/******/ \t\t// Check if module is in cache\n/******/ \t\tif(installedModules[moduleId])\n/******/ \t\t\treturn installedModules[moduleId].exports;\n\n/******/ \t\t// Create a new module (and put it into the cache)\n/******/ \t\tvar module = installedModules[moduleId] = {\n/******/ \t\t\texports: {},\n/******/ \t\t\tid: moduleId,\n/******/ \t\t\tloaded: false\n/******/ \t\t};\n\n/******/ \t\t// Execute the module function\n/******/ \t\tmodules[moduleId].call(module.exports, module, module.exports, __webpack_require__);\n\n/******/ \t\t// Flag the module as loaded\n/******/ \t\tmodule.loaded = true;\n\n/******/ \t\t// Return the exports of the module\n/******/ \t\treturn module.exports;\n/******/ \t}\n\n\n/******/ \t// expose the modules object (__webpack_modules__)\n/******/ \t__webpack_require__.m = modules;\n\n/******/ \t// expose the module cache\n/******/ \t__webpack_require__.c = installedModules;\n\n/******/ \t// __webpack_public_path__\n/******/ \t__webpack_require__.p = \"\";\n\n/******/ \t// Load entry module and return exports\n/******/ \treturn __webpack_require__(0);\n/******/ })\n/************************************************************************/\n/******/ ([\n/* 0 */\n/***/ function(module, exports, __webpack_require__) {\n\n\t//Use a partial d3 build, that doesn't need the DOM.\n\tvar d3f = __webpack_require__(1)\n\tvar Helpers = __webpack_require__(2)\n\n\n\tif (typeof onmessage !== 'undefined'){\n\t  onmessage = messaged\n\t} else {\n\t  function Wo(){}\n\t  Wo.prototype.postMessage = messaged\n\t  Wo.prototype.terminate = function(){}\n\t  module.exports = Wo\n\t}\n\n\tfunction messaged(event) {\n\t  var data = event\n\t  if (typeof event.data!== 'undefined') data = event.data\n\n\t  if (data.do === 'carto'){\n\n\t    var geo = data.geo,\n\t          topology = geo.topology,\n\t          geometries = geo.geometries,\n\t          path = geo.path,\n\t          translation = geo.translation,\n\t          projectionName = geo.projection.name,\n\t          scaling = geo.projection.scaling,\n\t          center = geo.projection.center,\n\t          translation = geo.projection.translation;\n\n\t    var values = data.values,\n\t        featureProperty = data.featureProperty,\n\t        task = data.task;\n\n\n\n\t    // copy it first\n\t    topology = Helpers.copy(topology);\n\n\n\n\t    // objects are projected into screen coordinates\n\t    // project the arcs into screen space\n\n\n\t    var projection =\n\t          d3f.geo[projectionName]()\n\t            .scale(scaling);\n\n\t    if (center != null) projection.center(center)\n\t    if (translation != null) projection.translate(translation)\n\n\t    var tf = Helpers.transformer(topology.transform),x,y,nArcVertices,vI,out1,nArcs=topology.arcs.length,aI=0,\n\t    projectedArcs = new Array(nArcs);\n\t    while(aI < nArcs){\n\t      x = 0;\n\t      y = 0;\n\t      nArcVertices = topology.arcs[aI].length;\n\t      vI = 0;\n\t      out1 = new Array(nArcVertices);\n\t      while( vI < nArcVertices){\n\t        topology.arcs[aI][vI][0] = (x += topology.arcs[aI][vI][0]);\n\t        topology.arcs[aI][vI][1] = (y += topology.arcs[aI][vI][1]);\n\t        out1[vI] = projection(tf(topology.arcs[aI][vI]));\n\t        vI++;\n\t      }\n\t      projectedArcs[aI++]=out1;\n\n\t    }\n\n\t    // path with identity projection\n\t    var path = d3f.geo.path()\n\t    .projection(null);\n\n\n\t    var objects = Helpers.object(projectedArcs, {type: \"GeometryCollection\", geometries: geometries})\n\t    .geometries.map(function(geom) {\n\t      return {\n\t        type: \"Feature\",\n\t        id: geom.id,\n\t        properties: Helpers.properties.call(null, geom, topology),\n\t        geometry: geom\n\t      };\n\t    });\n\n\t    function value(d){\n\t      return values[d.properties[featureProperty]]\n\t    }\n\n\t    var objectValues = objects.map(value),\n\t      totalValue = objectValues.reduce(function(a,b){return a + b;});\n\n\t    var iterations = 8;\n\n\n\t    //console.time(\"processing:\" + task)\n\t    var i = 0;\n\t    while (i++ < iterations) {\n\n\t      //var areas = objects.map(path.area)\n\t      //var totalArea = areas.reduce(function(a,b){return a + b}),\n\t      var areas = [], totalArea = 0;\n\t      for (var k = 0; k < objects.length; k++){\n\t        var area = path.area(objects[k])\n\t        areas.push(area)\n\t        totalArea += area\n\t      }\n\n\t      var sizeErrorsTot = 0,\n\t      sizeErrorsNum = 0;\n\n\t      ///for i = 1 to n do\n\t      var meta = []\n\t      for (var j = 0; j < objects.length; j++){\n\t        var o = objects[j],\n\t            area = Math.abs(areas[j]), // XXX: why do we have negative areas?\n\t            v = +objectValues[j],\n\t            ///Compute AD i , the desired area of the ith cell\n\t            desired = totalArea * v / totalValue,\n\t            radius = Math.sqrt(area / Math.PI),\n\t            mass = Math.sqrt(desired / Math.PI) - radius,\n\t            sizeError = Math.max(area, desired) / Math.min(area, desired);\n\n\t        sizeErrorsTot+=sizeError;\n\t        sizeErrorsNum++;\n\t        // console.log(o.id, \"@\", j, \"area:\", area, \"value:\", v, \"->\", desired, radius, mass, sizeError);\n\t        meta.push({\n\t          id:         o.id,\n\t          area:       area,\n\t          centroid:   path.centroid(o),\n\t          value:      v,\n\t          desired:    desired,\n\t          range: 100 * (Math.abs(desired - area)) / (Math.sqrt(Math.PI * area)),\n\t          radius:     radius,\n\t          mass:       mass,\n\t          sizeError:  sizeError\n\t        })\n\t      }\n\n\t      var sizeError = sizeErrorsTot / sizeErrorsNum,\n\t          forceReductionFactor = 1 / (1 + sizeError);\n\n\t      // console.log(\"meta:\", meta);\n\t      // console.log(\"  total area:\", totalArea);\n\t      // console.log(\"  force reduction factor:\", forceReductionFactor, \"mean error:\", sizeError);\n\n\t      var nArcVertices,vI,delta,nArcs=projectedArcs.length,aI=0,delta,nPolygon,pI,centroid,mass,radius,rSquared,dx,dy,distSquared,dist,Fij;\n\t      ///For each boundary line\n\t      while(aI < nArcs){\n\t        nArcVertices=projectedArcs[aI].length;\n\t        vI=0;\n\t        ///For each coordinate pair\n\t        while(vI < nArcVertices){\n\t          // create an array of vectors: [x, y]\n\t          delta = [0,0];\n\t          nPolygon = meta.length;\n\t          pI=0;\n\t          ///For each polygon centroid\n\t          while(pI < nPolygon) {\n\t            centroid =  meta[pI].centroid;\n\t            mass =      meta[pI].mass;\n\t            radius =    meta[pI].radius;\n\t            rSquared = (radius*radius);\n\t            dx = projectedArcs[aI][vI][0] - centroid[0];\n\t            dy = projectedArcs[aI][vI][1] - centroid[1];\n\t            distSquared = dx * dx + dy * dy;\n\t            dist=Math.sqrt(distSquared);\n\t            if (dist < meta[pI].range){\n\t              Fij = (dist > radius)\n\t              ? mass * radius / dist\n\t              : mass *\n\t              (distSquared / rSquared) *\n\t              (4 - 3 * dist / radius);\n\t              var tans = Helpers.arctans(dy, dx)\n\t              delta[0]+=(Fij * tans.cos);\n\t              delta[1]+=(Fij * tans.sin);\n\t            }\n\t            pI++;\n\t          }\n\t          projectedArcs[aI][vI][0] += (delta[0]*forceReductionFactor);\n\t          projectedArcs[aI][vI][1] += (delta[1]*forceReductionFactor);\n\t          vI++;\n\t        }\n\t        aI++;\n\t      }\n\n\t      // break if we hit the target size error\n\t      if (sizeError <= 1) break;\n\t    }\n\n\t    //console.timeEnd(\"processing:\" + task)\n\n\t    var response = {\n\t      done: 'processing',\n\t      //geohson featureCollection\n\t      features: objects,\n\t      //arcs can be useful to reconstruct topojson :\n\t      //arcs: projectedArcs,\n\t      task: task\n\t    }\n\n\t    if (typeof self !== 'undefined'){\n\t      self.postMessage(response)\n\t    } else {\n\t      this.onmessage(response)\n\t    }\n\t  }\n\t}\n\n\n/***/ },\n/* 1 */\n/***/ function(module, exports, __webpack_require__) {\n\n\tvar __WEBPACK_AMD_DEFINE_FACTORY__, __WEBPACK_AMD_DEFINE_RESULT__;!function(){\n\t  var d3f = {version: \"3.5.5\"}; // semver\n\tfunction d3f_identity(d) {\n\t  return d;\n\t}\n\tvar ε = 1e-6,\n\t    ε2 = ε * ε,\n\t    π = Math.PI,\n\t    τ = 2 * π,\n\t    τε = τ - ε,\n\t    halfπ = π / 2,\n\t    d3f_radians = π / 180,\n\t    d3f_degrees = 180 / π;\n\n\tfunction d3f_sgn(x) {\n\t  return x > 0 ? 1 : x < 0 ? -1 : 0;\n\t}\n\n\t// Returns the 2D cross product of AB and AC vectors, i.e., the z-component of\n\t// the 3D cross product in a quadrant I Cartesian coordinate system (+x is\n\t// right, +y is up). Returns a positive value if ABC is counter-clockwise,\n\t// negative if clockwise, and zero if the points are collinear.\n\tfunction d3f_cross2d(a, b, c) {\n\t  return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);\n\t}\n\n\tfunction d3f_acos(x) {\n\t  return x > 1 ? 0 : x < -1 ? π : Math.acos(x);\n\t}\n\n\tfunction d3f_asin(x) {\n\t  return x > 1 ? halfπ : x < -1 ? -halfπ : Math.asin(x);\n\t}\n\n\tfunction d3f_sinh(x) {\n\t  return ((x = Math.exp(x)) - 1 / x) / 2;\n\t}\n\n\tfunction d3f_cosh(x) {\n\t  return ((x = Math.exp(x)) + 1 / x) / 2;\n\t}\n\n\tfunction d3f_tanh(x) {\n\t  return ((x = Math.exp(2 * x)) - 1) / (x + 1);\n\t}\n\n\tfunction d3f_haversin(x) {\n\t  return (x = Math.sin(x / 2)) * x;\n\t}\n\td3f.geo = {};\n\t// Copies a variable number of methods from source to target.\n\td3f.rebind = function(target, source) {\n\t  var i = 1, n = arguments.length, method;\n\t  while (++i < n) target[method = arguments[i]] = d3f_rebind(target, source, source[method]);\n\t  return target;\n\t};\n\n\t// Method is assumed to be a standard D3 getter-setter:\n\t// If passed with no arguments, gets the value.\n\t// If passed with arguments, sets the value and returns the target.\n\tfunction d3f_rebind(target, source, method) {\n\t  return function() {\n\t    var value = method.apply(source, arguments);\n\t    return value === source ? target : value;\n\t  };\n\t}\n\tfunction d3f_true() {\n\t  return true;\n\t}\n\tvar abs = Math.abs;\n\td3f.merge = function(arrays) {\n\t  var n = arrays.length,\n\t      m,\n\t      i = -1,\n\t      j = 0,\n\t      merged,\n\t      array;\n\n\t  while (++i < n) j += arrays[i].length;\n\t  merged = new Array(j);\n\n\t  while (--n >= 0) {\n\t    array = arrays[n];\n\t    m = array.length;\n\t    while (--m >= 0) {\n\t      merged[--j] = array[m];\n\t    }\n\t  }\n\n\t  return merged;\n\t};\n\tfunction d3f_noop() {}\n\n\tfunction d3f_geo_spherical(cartesian) {\n\t  return [\n\t    Math.atan2(cartesian[1], cartesian[0]),\n\t    d3f_asin(cartesian[2])\n\t  ];\n\t}\n\n\tfunction d3f_geo_sphericalEqual(a, b) {\n\t  return abs(a[0] - b[0]) < ε && abs(a[1] - b[1]) < ε;\n\t}\n\n\t// General spherical polygon clipping algorithm: takes a polygon, cuts it into\n\t// visible line segments and rejoins the segments by interpolating along the\n\t// clip edge.\n\tfunction d3f_geo_clipPolygon(segments, compare, clipStartInside, interpolate, listener) {\n\t  var subject = [],\n\t      clip = [];\n\n\t  segments.forEach(function(segment) {\n\t    if ((n = segment.length - 1) <= 0) return;\n\t    var n, p0 = segment[0], p1 = segment[n];\n\n\t    // If the first and last points of a segment are coincident, then treat as\n\t    // a closed ring.\n\t    // TODO if all rings are closed, then the winding order of the exterior\n\t    // ring should be checked.\n\t    if (d3f_geo_sphericalEqual(p0, p1)) {\n\t      listener.lineStart();\n\t      for (var i = 0; i < n; ++i) listener.point((p0 = segment[i])[0], p0[1]);\n\t      listener.lineEnd();\n\t      return;\n\t    }\n\n\t    var a = new d3f_geo_clipPolygonIntersection(p0, segment, null, true),\n\t        b = new d3f_geo_clipPolygonIntersection(p0, null, a, false);\n\t    a.o = b;\n\t    subject.push(a);\n\t    clip.push(b);\n\t    a = new d3f_geo_clipPolygonIntersection(p1, segment, null, false);\n\t    b = new d3f_geo_clipPolygonIntersection(p1, null, a, true);\n\t    a.o = b;\n\t    subject.push(a);\n\t    clip.push(b);\n\t  });\n\t  clip.sort(compare);\n\t  d3f_geo_clipPolygonLinkCircular(subject);\n\t  d3f_geo_clipPolygonLinkCircular(clip);\n\t  if (!subject.length) return;\n\n\t  for (var i = 0, entry = clipStartInside, n = clip.length; i < n; ++i) {\n\t    clip[i].e = entry = !entry;\n\t  }\n\n\t  var start = subject[0],\n\t      points,\n\t      point;\n\t  while (1) {\n\t    // Find first unvisited intersection.\n\t    var current = start,\n\t        isSubject = true;\n\t    while (current.v) if ((current = current.n) === start) return;\n\t    points = current.z;\n\t    listener.lineStart();\n\t    do {\n\t      current.v = current.o.v = true;\n\t      if (current.e) {\n\t        if (isSubject) {\n\t          for (var i = 0, n = points.length; i < n; ++i) listener.point((point = points[i])[0], point[1]);\n\t        } else {\n\t          interpolate(current.x, current.n.x, 1, listener);\n\t        }\n\t        current = current.n;\n\t      } else {\n\t        if (isSubject) {\n\t          points = current.p.z;\n\t          for (var i = points.length - 1; i >= 0; --i) listener.point((point = points[i])[0], point[1]);\n\t        } else {\n\t          interpolate(current.x, current.p.x, -1, listener);\n\t        }\n\t        current = current.p;\n\t      }\n\t      current = current.o;\n\t      points = current.z;\n\t      isSubject = !isSubject;\n\t    } while (!current.v);\n\t    listener.lineEnd();\n\t  }\n\t}\n\n\tfunction d3f_geo_clipPolygonLinkCircular(array) {\n\t  if (!(n = array.length)) return;\n\t  var n,\n\t      i = 0,\n\t      a = array[0],\n\t      b;\n\t  while (++i < n) {\n\t    a.n = b = array[i];\n\t    b.p = a;\n\t    a = b;\n\t  }\n\t  a.n = b = array[0];\n\t  b.p = a;\n\t}\n\n\tfunction d3f_geo_clipPolygonIntersection(point, points, other, entry) {\n\t  this.x = point;\n\t  this.z = points;\n\t  this.o = other; // another intersection\n\t  this.e = entry; // is an entry?\n\t  this.v = false; // visited\n\t  this.n = this.p = null; // next & previous\n\t}\n\n\tfunction d3f_geo_clip(pointVisible, clipLine, interpolate, clipStart) {\n\t  return function(rotate, listener) {\n\t    var line = clipLine(listener),\n\t        rotatedClipStart = rotate.invert(clipStart[0], clipStart[1]);\n\n\t    var clip = {\n\t      point: point,\n\t      lineStart: lineStart,\n\t      lineEnd: lineEnd,\n\t      polygonStart: function() {\n\t        clip.point = pointRing;\n\t        clip.lineStart = ringStart;\n\t        clip.lineEnd = ringEnd;\n\t        segments = [];\n\t        polygon = [];\n\t      },\n\t      polygonEnd: function() {\n\t        clip.point = point;\n\t        clip.lineStart = lineStart;\n\t        clip.lineEnd = lineEnd;\n\n\t        segments = d3f.merge(segments);\n\t        var clipStartInside = d3f_geo_pointInPolygon(rotatedClipStart, polygon);\n\t        if (segments.length) {\n\t          if (!polygonStarted) listener.polygonStart(), polygonStarted = true;\n\t          d3f_geo_clipPolygon(segments, d3f_geo_clipSort, clipStartInside, interpolate, listener);\n\t        } else if (clipStartInside) {\n\t          if (!polygonStarted) listener.polygonStart(), polygonStarted = true;\n\t          listener.lineStart();\n\t          interpolate(null, null, 1, listener);\n\t          listener.lineEnd();\n\t        }\n\t        if (polygonStarted) listener.polygonEnd(), polygonStarted = false;\n\t        segments = polygon = null;\n\t      },\n\t      sphere: function() {\n\t        listener.polygonStart();\n\t        listener.lineStart();\n\t        interpolate(null, null, 1, listener);\n\t        listener.lineEnd();\n\t        listener.polygonEnd();\n\t      }\n\t    };\n\n\t    function point(λ, φ) {\n\t      var point = rotate(λ, φ);\n\t      if (pointVisible(λ = point[0], φ = point[1])) listener.point(λ, φ);\n\t    }\n\t    function pointLine(λ, φ) {\n\t      var point = rotate(λ, φ);\n\t      line.point(point[0], point[1]);\n\t    }\n\t    function lineStart() { clip.point = pointLine; line.lineStart(); }\n\t    function lineEnd() { clip.point = point; line.lineEnd(); }\n\n\t    var segments;\n\n\t    var buffer = d3f_geo_clipBufferListener(),\n\t        ringListener = clipLine(buffer),\n\t        polygonStarted = false,\n\t        polygon,\n\t        ring;\n\n\t    function pointRing(λ, φ) {\n\t      ring.push([λ, φ]);\n\t      var point = rotate(λ, φ);\n\t      ringListener.point(point[0], point[1]);\n\t    }\n\n\t    function ringStart() {\n\t      ringListener.lineStart();\n\t      ring = [];\n\t    }\n\n\t    function ringEnd() {\n\t      pointRing(ring[0][0], ring[0][1]);\n\t      ringListener.lineEnd();\n\n\t      var clean = ringListener.clean(),\n\t          ringSegments = buffer.buffer(),\n\t          segment,\n\t          n = ringSegments.length;\n\n\t      ring.pop();\n\t      polygon.push(ring);\n\t      ring = null;\n\n\t      if (!n) return;\n\n\t      // No intersections.\n\t      if (clean & 1) {\n\t        segment = ringSegments[0];\n\t        var n = segment.length - 1,\n\t            i = -1,\n\t            point;\n\t        if (n > 0) {\n\t          if (!polygonStarted) listener.polygonStart(), polygonStarted = true;\n\t          listener.lineStart();\n\t          while (++i < n) listener.point((point = segment[i])[0], point[1]);\n\t          listener.lineEnd();\n\t        }\n\t        return;\n\t      }\n\n\t      // Rejoin connected segments.\n\t      // TODO reuse bufferListener.rejoin()?\n\t      if (n > 1 && clean & 2) ringSegments.push(ringSegments.pop().concat(ringSegments.shift()));\n\n\t      segments.push(ringSegments.filter(d3f_geo_clipSegmentLength1));\n\t    }\n\n\t    return clip;\n\t  };\n\t}\n\n\tfunction d3f_geo_clipSegmentLength1(segment) {\n\t  return segment.length > 1;\n\t}\n\n\tfunction d3f_geo_clipBufferListener() {\n\t  var lines = [],\n\t      line;\n\t  return {\n\t    lineStart: function() { lines.push(line = []); },\n\t    point: function(λ, φ) { line.push([λ, φ]); },\n\t    lineEnd: d3f_noop,\n\t    buffer: function() {\n\t      var buffer = lines;\n\t      lines = [];\n\t      line = null;\n\t      return buffer;\n\t    },\n\t    rejoin: function() {\n\t      if (lines.length > 1) lines.push(lines.pop().concat(lines.shift()));\n\t    }\n\t  };\n\t}\n\n\t// Intersection points are sorted along the clip edge. For both antimeridian\n\t// cutting and circle clipping, the same comparison is used.\n\tfunction d3f_geo_clipSort(a, b) {\n\t  return ((a = a.x)[0] < 0 ? a[1] - halfπ - ε : halfπ - a[1])\n\t       - ((b = b.x)[0] < 0 ? b[1] - halfπ - ε : halfπ - b[1]);\n\t}\n\n\tvar d3f_geo_clipAntimeridian = d3f_geo_clip(\n\t    d3f_true,\n\t    d3f_geo_clipAntimeridianLine,\n\t    d3f_geo_clipAntimeridianInterpolate,\n\t    [-π, -π / 2]);\n\n\t// Takes a line and cuts into visible segments. Return values:\n\t//   0: there were intersections or the line was empty.\n\t//   1: no intersections.\n\t//   2: there were intersections, and the first and last segments should be\n\t//      rejoined.\n\tfunction d3f_geo_clipAntimeridianLine(listener) {\n\t  var λ0 = NaN,\n\t      φ0 = NaN,\n\t      sλ0 = NaN,\n\t      clean; // no intersections\n\n\t  return {\n\t    lineStart: function() {\n\t      listener.lineStart();\n\t      clean = 1;\n\t    },\n\t    point: function(λ1, φ1) {\n\t      var sλ1 = λ1 > 0 ? π : -π,\n\t          dλ = abs(λ1 - λ0);\n\t      if (abs(dλ - π) < ε) { // line crosses a pole\n\t        listener.point(λ0, φ0 = (φ0 + φ1) / 2 > 0 ? halfπ : -halfπ);\n\t        listener.point(sλ0, φ0);\n\t        listener.lineEnd();\n\t        listener.lineStart();\n\t        listener.point(sλ1, φ0);\n\t        listener.point(λ1, φ0);\n\t        clean = 0;\n\t      } else if (sλ0 !== sλ1 && dλ >= π) { // line crosses antimeridian\n\t        // handle degeneracies\n\t        if (abs(λ0 - sλ0) < ε) λ0 -= sλ0 * ε;\n\t        if (abs(λ1 - sλ1) < ε) λ1 -= sλ1 * ε;\n\t        φ0 = d3f_geo_clipAntimeridianIntersect(λ0, φ0, λ1, φ1);\n\t        listener.point(sλ0, φ0);\n\t        listener.lineEnd();\n\t        listener.lineStart();\n\t        listener.point(sλ1, φ0);\n\t        clean = 0;\n\t      }\n\t      listener.point(λ0 = λ1, φ0 = φ1);\n\t      sλ0 = sλ1;\n\t    },\n\t    lineEnd: function() {\n\t      listener.lineEnd();\n\t      λ0 = φ0 = NaN;\n\t    },\n\t    // if there are intersections, we always rejoin the first and last segments.\n\t    clean: function() { return 2 - clean; }\n\t  };\n\t}\n\n\tfunction d3f_geo_clipAntimeridianIntersect(λ0, φ0, λ1, φ1) {\n\t  var cosφ0,\n\t      cosφ1,\n\t      sinλ0_λ1 = Math.sin(λ0 - λ1);\n\t  return abs(sinλ0_λ1) > ε\n\t      ? Math.atan((Math.sin(φ0) * (cosφ1 = Math.cos(φ1)) * Math.sin(λ1)\n\t                 - Math.sin(φ1) * (cosφ0 = Math.cos(φ0)) * Math.sin(λ0))\n\t                 / (cosφ0 * cosφ1 * sinλ0_λ1))\n\t      : (φ0 + φ1) / 2;\n\t}\n\n\tfunction d3f_geo_clipAntimeridianInterpolate(from, to, direction, listener) {\n\t  var φ;\n\t  if (from == null) {\n\t    φ = direction * halfπ;\n\t    listener.point(-π,  φ);\n\t    listener.point( 0,  φ);\n\t    listener.point( π,  φ);\n\t    listener.point( π,  0);\n\t    listener.point( π, -φ);\n\t    listener.point( 0, -φ);\n\t    listener.point(-π, -φ);\n\t    listener.point(-π,  0);\n\t    listener.point(-π,  φ);\n\t  } else if (abs(from[0] - to[0]) > ε) {\n\t    var s = from[0] < to[0] ? π : -π;\n\t    φ = direction * s / 2;\n\t    listener.point(-s, φ);\n\t    listener.point( 0, φ);\n\t    listener.point( s, φ);\n\t  } else {\n\t    listener.point(to[0], to[1]);\n\t  }\n\t}\n\t// TODO\n\t// cross and scale return new vectors,\n\t// whereas add and normalize operate in-place\n\n\tfunction d3f_geo_cartesian(spherical) {\n\t  var λ = spherical[0],\n\t      φ = spherical[1],\n\t      cosφ = Math.cos(φ);\n\t  return [\n\t    cosφ * Math.cos(λ),\n\t    cosφ * Math.sin(λ),\n\t    Math.sin(φ)\n\t  ];\n\t}\n\n\tfunction d3f_geo_cartesianDot(a, b) {\n\t  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];\n\t}\n\n\tfunction d3f_geo_cartesianCross(a, b) {\n\t  return [\n\t    a[1] * b[2] - a[2] * b[1],\n\t    a[2] * b[0] - a[0] * b[2],\n\t    a[0] * b[1] - a[1] * b[0]\n\t  ];\n\t}\n\n\tfunction d3f_geo_cartesianAdd(a, b) {\n\t  a[0] += b[0];\n\t  a[1] += b[1];\n\t  a[2] += b[2];\n\t}\n\n\tfunction d3f_geo_cartesianScale(vector, k) {\n\t  return [\n\t    vector[0] * k,\n\t    vector[1] * k,\n\t    vector[2] * k\n\t  ];\n\t}\n\n\tfunction d3f_geo_cartesianNormalize(d) {\n\t  var l = Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);\n\t  d[0] /= l;\n\t  d[1] /= l;\n\t  d[2] /= l;\n\t}\n\tfunction d3f_geo_compose(a, b) {\n\n\t  function compose(x, y) {\n\t    return x = a(x, y), b(x[0], x[1]);\n\t  }\n\n\t  if (a.invert && b.invert) compose.invert = function(x, y) {\n\t    return x = b.invert(x, y), x && a.invert(x[0], x[1]);\n\t  };\n\n\t  return compose;\n\t}\n\n\tfunction d3f_geo_equirectangular(λ, φ) {\n\t  return [λ, φ];\n\t}\n\n\t(d3f.geo.equirectangular = function() {\n\t  return d3f_geo_projection(d3f_geo_equirectangular);\n\t}).raw = d3f_geo_equirectangular.invert = d3f_geo_equirectangular;\n\n\td3f.geo.rotation = function(rotate) {\n\t  rotate = d3f_geo_rotation(rotate[0] % 360 * d3f_radians, rotate[1] * d3f_radians, rotate.length > 2 ? rotate[2] * d3f_radians : 0);\n\n\t  function forward(coordinates) {\n\t    coordinates = rotate(coordinates[0] * d3f_radians, coordinates[1] * d3f_radians);\n\t    return coordinates[0] *= d3f_degrees, coordinates[1] *= d3f_degrees, coordinates;\n\t  }\n\n\t  forward.invert = function(coordinates) {\n\t    coordinates = rotate.invert(coordinates[0] * d3f_radians, coordinates[1] * d3f_radians);\n\t    return coordinates[0] *= d3f_degrees, coordinates[1] *= d3f_degrees, coordinates;\n\t  };\n\n\t  return forward;\n\t};\n\n\tfunction d3f_geo_identityRotation(λ, φ) {\n\t  return [λ > π ? λ - τ : λ < -π ? λ + τ : λ, φ];\n\t}\n\n\td3f_geo_identityRotation.invert = d3f_geo_equirectangular;\n\n\t// Note: |δλ| must be < 2π\n\tfunction d3f_geo_rotation(δλ, δφ, δγ) {\n\t  return δλ ? (δφ || δγ ? d3f_geo_compose(d3f_geo_rotationλ(δλ), d3f_geo_rotationφγ(δφ, δγ))\n\t    : d3f_geo_rotationλ(δλ))\n\t    : (δφ || δγ ? d3f_geo_rotationφγ(δφ, δγ)\n\t    : d3f_geo_identityRotation);\n\t}\n\n\tfunction d3f_geo_forwardRotationλ(δλ) {\n\t  return function(λ, φ) {\n\t    return λ += δλ, [λ > π ? λ - τ : λ < -π ? λ + τ : λ, φ];\n\t  };\n\t}\n\n\tfunction d3f_geo_rotationλ(δλ) {\n\t  var rotation = d3f_geo_forwardRotationλ(δλ);\n\t  rotation.invert = d3f_geo_forwardRotationλ(-δλ);\n\t  return rotation;\n\t}\n\n\tfunction d3f_geo_rotationφγ(δφ, δγ) {\n\t  var cosδφ = Math.cos(δφ),\n\t      sinδφ = Math.sin(δφ),\n\t      cosδγ = Math.cos(δγ),\n\t      sinδγ = Math.sin(δγ);\n\n\t  function rotation(λ, φ) {\n\t    var cosφ = Math.cos(φ),\n\t        x = Math.cos(λ) * cosφ,\n\t        y = Math.sin(λ) * cosφ,\n\t        z = Math.sin(φ),\n\t        k = z * cosδφ + x * sinδφ;\n\t    return [\n\t      Math.atan2(y * cosδγ - k * sinδγ, x * cosδφ - z * sinδφ),\n\t      d3f_asin(k * cosδγ + y * sinδγ)\n\t    ];\n\t  }\n\n\t  rotation.invert = function(λ, φ) {\n\t    var cosφ = Math.cos(φ),\n\t        x = Math.cos(λ) * cosφ,\n\t        y = Math.sin(λ) * cosφ,\n\t        z = Math.sin(φ),\n\t        k = z * cosδγ - y * sinδγ;\n\t    return [\n\t      Math.atan2(y * cosδγ + z * sinδγ, x * cosδφ + k * sinδφ),\n\t      d3f_asin(k * cosδφ - x * sinδφ)\n\t    ];\n\t  };\n\n\t  return rotation;\n\t}\n\n\td3f.geo.circle = function() {\n\t  var origin = [0, 0],\n\t      angle,\n\t      precision = 6,\n\t      interpolate;\n\n\t  function circle() {\n\t    var center = typeof origin === \"function\" ? origin.apply(this, arguments) : origin,\n\t        rotate = d3f_geo_rotation(-center[0] * d3f_radians, -center[1] * d3f_radians, 0).invert,\n\t        ring = [];\n\n\t    interpolate(null, null, 1, {\n\t      point: function(x, y) {\n\t        ring.push(x = rotate(x, y));\n\t        x[0] *= d3f_degrees, x[1] *= d3f_degrees;\n\t      }\n\t    });\n\n\t    return {type: \"Polygon\", coordinates: [ring]};\n\t  }\n\n\t  circle.origin = function(x) {\n\t    if (!arguments.length) return origin;\n\t    origin = x;\n\t    return circle;\n\t  };\n\n\t  circle.angle = function(x) {\n\t    if (!arguments.length) return angle;\n\t    interpolate = d3f_geo_circleInterpolate((angle = +x) * d3f_radians, precision * d3f_radians);\n\t    return circle;\n\t  };\n\n\t  circle.precision = function(_) {\n\t    if (!arguments.length) return precision;\n\t    interpolate = d3f_geo_circleInterpolate(angle * d3f_radians, (precision = +_) * d3f_radians);\n\t    return circle;\n\t  };\n\n\t  return circle.angle(90);\n\t};\n\n\t// Interpolates along a circle centered at [0°, 0°], with a given radius and\n\t// precision.\n\tfunction d3f_geo_circleInterpolate(radius, precision) {\n\t  var cr = Math.cos(radius),\n\t      sr = Math.sin(radius);\n\t  return function(from, to, direction, listener) {\n\t    var step = direction * precision;\n\t    if (from != null) {\n\t      from = d3f_geo_circleAngle(cr, from);\n\t      to = d3f_geo_circleAngle(cr, to);\n\t      if (direction > 0 ? from < to: from > to) from += direction * τ;\n\t    } else {\n\t      from = radius + direction * τ;\n\t      to = radius - .5 * step;\n\t    }\n\t    for (var point, t = from; direction > 0 ? t > to : t < to; t -= step) {\n\t      listener.point((point = d3f_geo_spherical([\n\t        cr,\n\t        -sr * Math.cos(t),\n\t        -sr * Math.sin(t)\n\t      ]))[0], point[1]);\n\t    }\n\t  };\n\t}\n\n\t// Signed angle of a cartesian point relative to [cr, 0, 0].\n\tfunction d3f_geo_circleAngle(cr, point) {\n\t  var a = d3f_geo_cartesian(point);\n\t  a[0] -= cr;\n\t  d3f_geo_cartesianNormalize(a);\n\t  var angle = d3f_acos(-a[1]);\n\t  return ((-a[2] < 0 ? -angle : angle) + 2 * Math.PI - ε) % (2 * Math.PI);\n\t}\n\t// Adds floating point numbers with twice the normal precision.\n\t// Reference: J. R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and\n\t// Fast Robust Geometric Predicates, Discrete & Computational Geometry 18(3)\n\t// 305–363 (1997).\n\t// Code adapted from GeographicLib by Charles F. F. Karney,\n\t// http://geographiclib.sourceforge.net/\n\t// See lib/geographiclib/LICENSE for details.\n\n\tfunction d3f_adder() {}\n\n\td3f_adder.prototype = {\n\t  s: 0, // rounded value\n\t  t: 0, // exact error\n\t  add: function(y) {\n\t    d3f_adderSum(y, this.t, d3f_adderTemp);\n\t    d3f_adderSum(d3f_adderTemp.s, this.s, this);\n\t    if (this.s) this.t += d3f_adderTemp.t;\n\t    else this.s = d3f_adderTemp.t;\n\t  },\n\t  reset: function() {\n\t    this.s = this.t = 0;\n\t  },\n\t  valueOf: function() {\n\t    return this.s;\n\t  }\n\t};\n\n\tvar d3f_adderTemp = new d3f_adder;\n\n\tfunction d3f_adderSum(a, b, o) {\n\t  var x = o.s = a + b, // a + b\n\t      bv = x - a, av = x - bv; // b_virtual & a_virtual\n\t  o.t = (a - av) + (b - bv); // a_roundoff + b_roundoff\n\t}\n\n\td3f.geo.stream = function(object, listener) {\n\t  if (object && d3f_geo_streamObjectType.hasOwnProperty(object.type)) {\n\t    d3f_geo_streamObjectType[object.type](object, listener);\n\t  } else {\n\t    d3f_geo_streamGeometry(object, listener);\n\t  }\n\t};\n\n\tfunction d3f_geo_streamGeometry(geometry, listener) {\n\t  if (geometry && d3f_geo_streamGeometryType.hasOwnProperty(geometry.type)) {\n\t    d3f_geo_streamGeometryType[geometry.type](geometry, listener);\n\t  }\n\t}\n\n\tvar d3f_geo_streamObjectType = {\n\t  Feature: function(feature, listener) {\n\t    d3f_geo_streamGeometry(feature.geometry, listener);\n\t  },\n\t  FeatureCollection: function(object, listener) {\n\t    var features = object.features, i = -1, n = features.length;\n\t    while (++i < n) d3f_geo_streamGeometry(features[i].geometry, listener);\n\t  }\n\t};\n\n\tvar d3f_geo_streamGeometryType = {\n\t  Sphere: function(object, listener) {\n\t    listener.sphere();\n\t  },\n\t  Point: function(object, listener) {\n\t    object = object.coordinates;\n\t    listener.point(object[0], object[1], object[2]);\n\t  },\n\t  MultiPoint: function(object, listener) {\n\t    var coordinates = object.coordinates, i = -1, n = coordinates.length;\n\t    while (++i < n) object = coordinates[i], listener.point(object[0], object[1], object[2]);\n\t  },\n\t  LineString: function(object, listener) {\n\t    d3f_geo_streamLine(object.coordinates, listener, 0);\n\t  },\n\t  MultiLineString: function(object, listener) {\n\t    var coordinates = object.coordinates, i = -1, n = coordinates.length;\n\t    while (++i < n) d3f_geo_streamLine(coordinates[i], listener, 0);\n\t  },\n\t  Polygon: function(object, listener) {\n\t    d3f_geo_streamPolygon(object.coordinates, listener);\n\t  },\n\t  MultiPolygon: function(object, listener) {\n\t    var coordinates = object.coordinates, i = -1, n = coordinates.length;\n\t    while (++i < n) d3f_geo_streamPolygon(coordinates[i], listener);\n\t  },\n\t  GeometryCollection: function(object, listener) {\n\t    var geometries = object.geometries, i = -1, n = geometries.length;\n\t    while (++i < n) d3f_geo_streamGeometry(geometries[i], listener);\n\t  }\n\t};\n\n\tfunction d3f_geo_streamLine(coordinates, listener, closed) {\n\t  var i = -1, n = coordinates.length - closed, coordinate;\n\t  listener.lineStart();\n\t  while (++i < n) coordinate = coordinates[i], listener.point(coordinate[0], coordinate[1], coordinate[2]);\n\t  listener.lineEnd();\n\t}\n\n\tfunction d3f_geo_streamPolygon(coordinates, listener) {\n\t  var i = -1, n = coordinates.length;\n\t  listener.polygonStart();\n\t  while (++i < n) d3f_geo_streamLine(coordinates[i], listener, 1);\n\t  listener.polygonEnd();\n\t}\n\n\td3f.geo.area = function(object) {\n\t  d3f_geo_areaSum = 0;\n\t  d3f.geo.stream(object, d3f_geo_area);\n\t  return d3f_geo_areaSum;\n\t};\n\n\tvar d3f_geo_areaSum,\n\t    d3f_geo_areaRingSum = new d3f_adder;\n\n\tvar d3f_geo_area = {\n\t  sphere: function() { d3f_geo_areaSum += 4 * π; },\n\t  point: d3f_noop,\n\t  lineStart: d3f_noop,\n\t  lineEnd: d3f_noop,\n\n\t  // Only count area for polygon rings.\n\t  polygonStart: function() {\n\t    d3f_geo_areaRingSum.reset();\n\t    d3f_geo_area.lineStart = d3f_geo_areaRingStart;\n\t  },\n\t  polygonEnd: function() {\n\t    var area = 2 * d3f_geo_areaRingSum;\n\t    d3f_geo_areaSum += area < 0 ? 4 * π + area : area;\n\t    d3f_geo_area.lineStart = d3f_geo_area.lineEnd = d3f_geo_area.point = d3f_noop;\n\t  }\n\t};\n\n\tfunction d3f_geo_areaRingStart() {\n\t  var λ00, φ00, λ0, cosφ0, sinφ0; // start point and previous point\n\n\t  // For the first point, …\n\t  d3f_geo_area.point = function(λ, φ) {\n\t    d3f_geo_area.point = nextPoint;\n\t    λ0 = (λ00 = λ) * d3f_radians, cosφ0 = Math.cos(φ = (φ00 = φ) * d3f_radians / 2 + π / 4), sinφ0 = Math.sin(φ);\n\t  };\n\n\t  // For subsequent points, …\n\t  function nextPoint(λ, φ) {\n\t    λ *= d3f_radians;\n\t    φ = φ * d3f_radians / 2 + π / 4; // half the angular distance from south pole\n\n\t    // Spherical excess E for a spherical triangle with vertices: south pole,\n\t    // previous point, current point.  Uses a formula derived from Cagnoli’s\n\t    // theorem.  See Todhunter, Spherical Trig. (1871), Sec. 103, Eq. (2).\n\t    var dλ = λ - λ0,\n\t        sdλ = dλ >= 0 ? 1 : -1,\n\t        adλ = sdλ * dλ,\n\t        cosφ = Math.cos(φ),\n\t        sinφ = Math.sin(φ),\n\t        k = sinφ0 * sinφ,\n\t        u = cosφ0 * cosφ + k * Math.cos(adλ),\n\t        v = k * sdλ * Math.sin(adλ);\n\t    d3f_geo_areaRingSum.add(Math.atan2(v, u));\n\n\t    // Advance the previous points.\n\t    λ0 = λ, cosφ0 = cosφ, sinφ0 = sinφ;\n\t  }\n\n\t  // For the last point, return to the start.\n\t  d3f_geo_area.lineEnd = function() {\n\t    nextPoint(λ00, φ00);\n\t  };\n\t}\n\n\tfunction d3f_geo_pointInPolygon(point, polygon) {\n\t  var meridian = point[0],\n\t      parallel = point[1],\n\t      meridianNormal = [Math.sin(meridian), -Math.cos(meridian), 0],\n\t      polarAngle = 0,\n\t      winding = 0;\n\t  d3f_geo_areaRingSum.reset();\n\n\t  for (var i = 0, n = polygon.length; i < n; ++i) {\n\t    var ring = polygon[i],\n\t        m = ring.length;\n\t    if (!m) continue;\n\t    var point0 = ring[0],\n\t        λ0 = point0[0],\n\t        φ0 = point0[1] / 2 + π / 4,\n\t        sinφ0 = Math.sin(φ0),\n\t        cosφ0 = Math.cos(φ0),\n\t        j = 1;\n\n\t    while (true) {\n\t      if (j === m) j = 0;\n\t      point = ring[j];\n\t      var λ = point[0],\n\t          φ = point[1] / 2 + π / 4,\n\t          sinφ = Math.sin(φ),\n\t          cosφ = Math.cos(φ),\n\t          dλ = λ - λ0,\n\t          sdλ = dλ >= 0 ? 1 : -1,\n\t          adλ = sdλ * dλ,\n\t          antimeridian = adλ > π,\n\t          k = sinφ0 * sinφ;\n\t      d3f_geo_areaRingSum.add(Math.atan2(k * sdλ * Math.sin(adλ), cosφ0 * cosφ + k * Math.cos(adλ)));\n\n\t      polarAngle += antimeridian ? dλ + sdλ * τ : dλ;\n\n\t      // Are the longitudes either side of the point's meridian, and are the\n\t      // latitudes smaller than the parallel?\n\t      if (antimeridian ^ λ0 >= meridian ^ λ >= meridian) {\n\t        var arc = d3f_geo_cartesianCross(d3f_geo_cartesian(point0), d3f_geo_cartesian(point));\n\t        d3f_geo_cartesianNormalize(arc);\n\t        var intersection = d3f_geo_cartesianCross(meridianNormal, arc);\n\t        d3f_geo_cartesianNormalize(intersection);\n\t        var φarc = (antimeridian ^ dλ >= 0 ? -1 : 1) * d3f_asin(intersection[2]);\n\t        if (parallel > φarc || parallel === φarc && (arc[0] || arc[1])) {\n\t          winding += antimeridian ^ dλ >= 0 ? 1 : -1;\n\t        }\n\t      }\n\t      if (!j++) break;\n\t      λ0 = λ, sinφ0 = sinφ, cosφ0 = cosφ, point0 = point;\n\t    }\n\t  }\n\n\t  // First, determine whether the South pole is inside or outside:\n\t  //\n\t  // It is inside if:\n\t  // * the polygon winds around it in a clockwise direction.\n\t  // * the polygon does not (cumulatively) wind around it, but has a negative\n\t  //   (counter-clockwise) area.\n\t  //\n\t  // Second, count the (signed) number of times a segment crosses a meridian\n\t  // from the point to the South pole.  If it is zero, then the point is the\n\t  // same side as the South pole.\n\n\t  return (polarAngle < -ε || polarAngle < ε && d3f_geo_areaRingSum < 0) ^ (winding & 1);\n\t}\n\n\t// Clip features against a small circle centered at [0°, 0°].\n\tfunction d3f_geo_clipCircle(radius) {\n\t  var cr = Math.cos(radius),\n\t      smallRadius = cr > 0,\n\t      notHemisphere = abs(cr) > ε, // TODO optimise for this common case\n\t      interpolate = d3f_geo_circleInterpolate(radius, 6 * d3f_radians);\n\n\t  return d3f_geo_clip(visible, clipLine, interpolate, smallRadius ? [0, -radius] : [-π, radius - π]);\n\n\t  function visible(λ, φ) {\n\t    return Math.cos(λ) * Math.cos(φ) > cr;\n\t  }\n\n\t  // Takes a line and cuts into visible segments. Return values used for\n\t  // polygon clipping:\n\t  //   0: there were intersections or the line was empty.\n\t  //   1: no intersections.\n\t  //   2: there were intersections, and the first and last segments should be\n\t  //      rejoined.\n\t  function clipLine(listener) {\n\t    var point0, // previous point\n\t        c0, // code for previous point\n\t        v0, // visibility of previous point\n\t        v00, // visibility of first point\n\t        clean; // no intersections\n\t    return {\n\t      lineStart: function() {\n\t        v00 = v0 = false;\n\t        clean = 1;\n\t      },\n\t      point: function(λ, φ) {\n\t        var point1 = [λ, φ],\n\t            point2,\n\t            v = visible(λ, φ),\n\t            c = smallRadius\n\t              ? v ? 0 : code(λ, φ)\n\t              : v ? code(λ + (λ < 0 ? π : -π), φ) : 0;\n\t        if (!point0 && (v00 = v0 = v)) listener.lineStart();\n\t        // Handle degeneracies.\n\t        // TODO ignore if not clipping polygons.\n\t        if (v !== v0) {\n\t          point2 = intersect(point0, point1);\n\t          if (d3f_geo_sphericalEqual(point0, point2) || d3f_geo_sphericalEqual(point1, point2)) {\n\t            point1[0] += ε;\n\t            point1[1] += ε;\n\t            v = visible(point1[0], point1[1]);\n\t          }\n\t        }\n\t        if (v !== v0) {\n\t          clean = 0;\n\t          if (v) {\n\t            // outside going in\n\t            listener.lineStart();\n\t            point2 = intersect(point1, point0);\n\t            listener.point(point2[0], point2[1]);\n\t          } else {\n\t            // inside going out\n\t            point2 = intersect(point0, point1);\n\t            listener.point(point2[0], point2[1]);\n\t            listener.lineEnd();\n\t          }\n\t          point0 = point2;\n\t        } else if (notHemisphere && point0 && smallRadius ^ v) {\n\t          var t;\n\t          // If the codes for two points are different, or are both zero,\n\t          // and there this segment intersects with the small circle.\n\t          if (!(c & c0) && (t = intersect(point1, point0, true))) {\n\t            clean = 0;\n\t            if (smallRadius) {\n\t              listener.lineStart();\n\t              listener.point(t[0][0], t[0][1]);\n\t              listener.point(t[1][0], t[1][1]);\n\t              listener.lineEnd();\n\t            } else {\n\t              listener.point(t[1][0], t[1][1]);\n\t              listener.lineEnd();\n\t              listener.lineStart();\n\t              listener.point(t[0][0], t[0][1]);\n\t            }\n\t          }\n\t        }\n\t        if (v && (!point0 || !d3f_geo_sphericalEqual(point0, point1))) {\n\t          listener.point(point1[0], point1[1]);\n\t        }\n\t        point0 = point1, v0 = v, c0 = c;\n\t      },\n\t      lineEnd: function() {\n\t        if (v0) listener.lineEnd();\n\t        point0 = null;\n\t      },\n\t      // Rejoin first and last segments if there were intersections and the first\n\t      // and last points were visible.\n\t      clean: function() { return clean | ((v00 && v0) << 1); }\n\t    };\n\t  }\n\n\t  // Intersects the great circle between a and b with the clip circle.\n\t  function intersect(a, b, two) {\n\t    var pa = d3f_geo_cartesian(a),\n\t        pb = d3f_geo_cartesian(b);\n\n\t    // We have two planes, n1.p = d1 and n2.p = d2.\n\t    // Find intersection line p(t) = c1 n1 + c2 n2 + t (n1 ⨯ n2).\n\t    var n1 = [1, 0, 0], // normal\n\t        n2 = d3f_geo_cartesianCross(pa, pb),\n\t        n2n2 = d3f_geo_cartesianDot(n2, n2),\n\t        n1n2 = n2[0], // d3f_geo_cartesianDot(n1, n2),\n\t        determinant = n2n2 - n1n2 * n1n2;\n\n\t    // Two polar points.\n\t    if (!determinant) return !two && a;\n\n\t    var c1 =  cr * n2n2 / determinant,\n\t        c2 = -cr * n1n2 / determinant,\n\t        n1xn2 = d3f_geo_cartesianCross(n1, n2),\n\t        A = d3f_geo_cartesianScale(n1, c1),\n\t        B = d3f_geo_cartesianScale(n2, c2);\n\t    d3f_geo_cartesianAdd(A, B);\n\n\t    // Solve |p(t)|^2 = 1.\n\t    var u = n1xn2,\n\t        w = d3f_geo_cartesianDot(A, u),\n\t        uu = d3f_geo_cartesianDot(u, u),\n\t        t2 = w * w - uu * (d3f_geo_cartesianDot(A, A) - 1);\n\n\t    if (t2 < 0) return;\n\n\t    var t = Math.sqrt(t2),\n\t        q = d3f_geo_cartesianScale(u, (-w - t) / uu);\n\t    d3f_geo_cartesianAdd(q, A);\n\t    q = d3f_geo_spherical(q);\n\t    if (!two) return q;\n\n\t    // Two intersection points.\n\t    var λ0 = a[0],\n\t        λ1 = b[0],\n\t        φ0 = a[1],\n\t        φ1 = b[1],\n\t        z;\n\t    if (λ1 < λ0) z = λ0, λ0 = λ1, λ1 = z;\n\t    var δλ = λ1 - λ0,\n\t        polar = abs(δλ - π) < ε,\n\t        meridian = polar || δλ < ε;\n\n\t    if (!polar && φ1 < φ0) z = φ0, φ0 = φ1, φ1 = z;\n\n\t    // Check that the first point is between a and b.\n\t    if (meridian\n\t        ? polar\n\t          ? φ0 + φ1 > 0 ^ q[1] < (abs(q[0] - λ0) < ε ? φ0 : φ1)\n\t          : φ0 <= q[1] && q[1] <= φ1\n\t        : δλ > π ^ (λ0 <= q[0] && q[0] <= λ1)) {\n\t      var q1 = d3f_geo_cartesianScale(u, (-w + t) / uu);\n\t      d3f_geo_cartesianAdd(q1, A);\n\t      return [q, d3f_geo_spherical(q1)];\n\t    }\n\t  }\n\n\t  // Generates a 4-bit vector representing the location of a point relative to\n\t  // the small circle's bounding box.\n\t  function code(λ, φ) {\n\t    var r = smallRadius ? radius : π - radius,\n\t        code = 0;\n\t    if (λ < -r) code |= 1; // left\n\t    else if (λ > r) code |= 2; // right\n\t    if (φ < -r) code |= 4; // below\n\t    else if (φ > r) code |= 8; // above\n\t    return code;\n\t  }\n\t}\n\n\t// Liang–Barsky line clipping.\n\tfunction d3f_geom_clipLine(x0, y0, x1, y1) {\n\t  return function(line) {\n\t    var a = line.a,\n\t        b = line.b,\n\t        ax = a.x,\n\t        ay = a.y,\n\t        bx = b.x,\n\t        by = b.y,\n\t        t0 = 0,\n\t        t1 = 1,\n\t        dx = bx - ax,\n\t        dy = by - ay,\n\t        r;\n\n\t    r = x0 - ax;\n\t    if (!dx && r > 0) return;\n\t    r /= dx;\n\t    if (dx < 0) {\n\t      if (r < t0) return;\n\t      if (r < t1) t1 = r;\n\t    } else if (dx > 0) {\n\t      if (r > t1) return;\n\t      if (r > t0) t0 = r;\n\t    }\n\n\t    r = x1 - ax;\n\t    if (!dx && r < 0) return;\n\t    r /= dx;\n\t    if (dx < 0) {\n\t      if (r > t1) return;\n\t      if (r > t0) t0 = r;\n\t    } else if (dx > 0) {\n\t      if (r < t0) return;\n\t      if (r < t1) t1 = r;\n\t    }\n\n\t    r = y0 - ay;\n\t    if (!dy && r > 0) return;\n\t    r /= dy;\n\t    if (dy < 0) {\n\t      if (r < t0) return;\n\t      if (r < t1) t1 = r;\n\t    } else if (dy > 0) {\n\t      if (r > t1) return;\n\t      if (r > t0) t0 = r;\n\t    }\n\n\t    r = y1 - ay;\n\t    if (!dy && r < 0) return;\n\t    r /= dy;\n\t    if (dy < 0) {\n\t      if (r > t1) return;\n\t      if (r > t0) t0 = r;\n\t    } else if (dy > 0) {\n\t      if (r < t0) return;\n\t      if (r < t1) t1 = r;\n\t    }\n\n\t    if (t0 > 0) line.a = {x: ax + t0 * dx, y: ay + t0 * dy};\n\t    if (t1 < 1) line.b = {x: ax + t1 * dx, y: ay + t1 * dy};\n\t    return line;\n\t  };\n\t}\n\n\tvar d3f_geo_clipExtentMAX = 1e9;\n\n\td3f.geo.clipExtent = function() {\n\t  var x0, y0, x1, y1,\n\t      stream,\n\t      clip,\n\t      clipExtent = {\n\t        stream: function(output) {\n\t          if (stream) stream.valid = false;\n\t          stream = clip(output);\n\t          stream.valid = true; // allow caching by d3f.geo.path\n\t          return stream;\n\t        },\n\t        extent: function(_) {\n\t          if (!arguments.length) return [[x0, y0], [x1, y1]];\n\t          clip = d3f_geo_clipExtent(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]);\n\t          if (stream) stream.valid = false, stream = null;\n\t          return clipExtent;\n\t        }\n\t      };\n\t  return clipExtent.extent([[0, 0], [960, 500]]);\n\t};\n\n\tfunction d3f_geo_clipExtent(x0, y0, x1, y1) {\n\t  return function(listener) {\n\t    var listener_ = listener,\n\t        bufferListener = d3f_geo_clipBufferListener(),\n\t        clipLine = d3f_geom_clipLine(x0, y0, x1, y1),\n\t        segments,\n\t        polygon,\n\t        ring;\n\n\t    var clip = {\n\t      point: point,\n\t      lineStart: lineStart,\n\t      lineEnd: lineEnd,\n\t      polygonStart: function() {\n\t        listener = bufferListener;\n\t        segments = [];\n\t        polygon = [];\n\t        clean = true;\n\t      },\n\t      polygonEnd: function() {\n\t        listener = listener_;\n\t        segments = d3f.merge(segments);\n\t        var clipStartInside = insidePolygon([x0, y1]),\n\t            inside = clean && clipStartInside,\n\t            visible = segments.length;\n\t        if (inside || visible) {\n\t          listener.polygonStart();\n\t          if (inside) {\n\t            listener.lineStart();\n\t            interpolate(null, null, 1, listener);\n\t            listener.lineEnd();\n\t          }\n\t          if (visible) {\n\t            d3f_geo_clipPolygon(segments, compare, clipStartInside, interpolate, listener);\n\t          }\n\t          listener.polygonEnd();\n\t        }\n\t        segments = polygon = ring = null;\n\t      }\n\t    };\n\n\t    function insidePolygon(p) {\n\t      var wn = 0, // the winding number counter\n\t          n = polygon.length,\n\t          y = p[1];\n\n\t      for (var i = 0; i < n; ++i) {\n\t        for (var j = 1, v = polygon[i], m = v.length, a = v[0], b; j < m; ++j) {\n\t          b = v[j];\n\t          if (a[1] <= y) {\n\t            if (b[1] >  y && d3f_cross2d(a, b, p) > 0) ++wn;\n\t          } else {\n\t            if (b[1] <= y && d3f_cross2d(a, b, p) < 0) --wn;\n\t          }\n\t          a = b;\n\t        }\n\t      }\n\t      return wn !== 0;\n\t    }\n\n\t    function interpolate(from, to, direction, listener) {\n\t      var a = 0, a1 = 0;\n\t      if (from == null ||\n\t          (a = corner(from, direction)) !== (a1 = corner(to, direction)) ||\n\t          comparePoints(from, to) < 0 ^ direction > 0) {\n\t        do {\n\t          listener.point(a === 0 || a === 3 ? x0 : x1, a > 1 ? y1 : y0);\n\t        } while ((a = (a + direction + 4) % 4) !== a1);\n\t      } else {\n\t        listener.point(to[0], to[1]);\n\t      }\n\t    }\n\n\t    function pointVisible(x, y) {\n\t      return x0 <= x && x <= x1 && y0 <= y && y <= y1;\n\t    }\n\n\t    function point(x, y) {\n\t      if (pointVisible(x, y)) listener.point(x, y);\n\t    }\n\n\t    var x__, y__, v__, // first point\n\t        x_, y_, v_, // previous point\n\t        first,\n\t        clean;\n\n\t    function lineStart() {\n\t      clip.point = linePoint;\n\t      if (polygon) polygon.push(ring = []);\n\t      first = true;\n\t      v_ = false;\n\t      x_ = y_ = NaN;\n\t    }\n\n\t    function lineEnd() {\n\t      // TODO rather than special-case polygons, simply handle them separately.\n\t      // Ideally, coincident intersection points should be jittered to avoid\n\t      // clipping issues.\n\t      if (segments) {\n\t        linePoint(x__, y__);\n\t        if (v__ && v_) bufferListener.rejoin();\n\t        segments.push(bufferListener.buffer());\n\t      }\n\t      clip.point = point;\n\t      if (v_) listener.lineEnd();\n\t    }\n\n\t    function linePoint(x, y) {\n\t      x = Math.max(-d3f_geo_clipExtentMAX, Math.min(d3f_geo_clipExtentMAX, x));\n\t      y = Math.max(-d3f_geo_clipExtentMAX, Math.min(d3f_geo_clipExtentMAX, y));\n\t      var v = pointVisible(x, y);\n\t      if (polygon) ring.push([x, y]);\n\t      if (first) {\n\t        x__ = x, y__ = y, v__ = v;\n\t        first = false;\n\t        if (v) {\n\t          listener.lineStart();\n\t          listener.point(x, y);\n\t        }\n\t      } else {\n\t        if (v && v_) listener.point(x, y);\n\t        else {\n\t          var l = {a: {x: x_, y: y_}, b: {x: x, y: y}};\n\t          if (clipLine(l)) {\n\t            if (!v_) {\n\t              listener.lineStart();\n\t              listener.point(l.a.x, l.a.y);\n\t            }\n\t            listener.point(l.b.x, l.b.y);\n\t            if (!v) listener.lineEnd();\n\t            clean = false;\n\t          } else if (v) {\n\t            listener.lineStart();\n\t            listener.point(x, y);\n\t            clean = false;\n\t          }\n\t        }\n\t      }\n\t      x_ = x, y_ = y, v_ = v;\n\t    }\n\n\t    return clip;\n\t  };\n\n\t  function corner(p, direction) {\n\t    return abs(p[0] - x0) < ε ? direction > 0 ? 0 : 3\n\t        : abs(p[0] - x1) < ε ? direction > 0 ? 2 : 1\n\t        : abs(p[1] - y0) < ε ? direction > 0 ? 1 : 0\n\t        : direction > 0 ? 3 : 2; // abs(p[1] - y1) < ε\n\t  }\n\n\t  function compare(a, b) {\n\t    return comparePoints(a.x, b.x);\n\t  }\n\n\t  function comparePoints(a, b) {\n\t    var ca = corner(a, 1),\n\t        cb = corner(b, 1);\n\t    return ca !== cb ? ca - cb\n\t        : ca === 0 ? b[1] - a[1]\n\t        : ca === 1 ? a[0] - b[0]\n\t        : ca === 2 ? a[1] - b[1]\n\t        : b[0] - a[0];\n\t  }\n\t}\n\n\tfunction d3f_geo_resample(project) {\n\t  var δ2 = .5, // precision, px²\n\t      cosMinDistance = Math.cos(30 * d3f_radians), // cos(minimum angular distance)\n\t      maxDepth = 16;\n\n\t  function resample(stream) {\n\t    return (maxDepth ? resampleRecursive : resampleNone)(stream);\n\t  }\n\n\t  function resampleNone(stream) {\n\t    return d3f_geo_transformPoint(stream, function(x, y) {\n\t      x = project(x, y);\n\t      stream.point(x[0], x[1]);\n\t    });\n\t  }\n\n\t  function resampleRecursive(stream) {\n\t    var λ00, φ00, x00, y00, a00, b00, c00, // first point\n\t        λ0, x0, y0, a0, b0, c0; // previous point\n\n\t    var resample = {\n\t      point: point,\n\t      lineStart: lineStart,\n\t      lineEnd: lineEnd,\n\t      polygonStart: function() { stream.polygonStart(); resample.lineStart = ringStart; },\n\t      polygonEnd: function() { stream.polygonEnd(); resample.lineStart = lineStart; }\n\t    };\n\n\t    function point(x, y) {\n\t      x = project(x, y);\n\t      stream.point(x[0], x[1]);\n\t    }\n\n\t    function lineStart() {\n\t      x0 = NaN;\n\t      resample.point = linePoint;\n\t      stream.lineStart();\n\t    }\n\n\t    function linePoint(λ, φ) {\n\t      var c = d3f_geo_cartesian([λ, φ]), p = project(λ, φ);\n\t      resampleLineTo(x0, y0, λ0, a0, b0, c0, x0 = p[0], y0 = p[1], λ0 = λ, a0 = c[0], b0 = c[1], c0 = c[2], maxDepth, stream);\n\t      stream.point(x0, y0);\n\t    }\n\n\t    function lineEnd() {\n\t      resample.point = point;\n\t      stream.lineEnd();\n\t    }\n\n\t    function ringStart() {\n\t      lineStart();\n\t      resample.point = ringPoint;\n\t      resample.lineEnd = ringEnd;\n\t    }\n\n\t    function ringPoint(λ, φ) {\n\t      linePoint(λ00 = λ, φ00 = φ), x00 = x0, y00 = y0, a00 = a0, b00 = b0, c00 = c0;\n\t      resample.point = linePoint;\n\t    }\n\n\t    function ringEnd() {\n\t      resampleLineTo(x0, y0, λ0, a0, b0, c0, x00, y00, λ00, a00, b00, c00, maxDepth, stream);\n\t      resample.lineEnd = lineEnd;\n\t      lineEnd();\n\t    }\n\n\t    return resample;\n\t  }\n\n\t  function resampleLineTo(x0, y0, λ0, a0, b0, c0, x1, y1, λ1, a1, b1, c1, depth, stream) {\n\t    var dx = x1 - x0,\n\t        dy = y1 - y0,\n\t        d2 = dx * dx + dy * dy;\n\t    if (d2 > 4 * δ2 && depth--) {\n\t      var a = a0 + a1,\n\t          b = b0 + b1,\n\t          c = c0 + c1,\n\t          m = Math.sqrt(a * a + b * b + c * c),\n\t          φ2 = Math.asin(c /= m),\n\t          λ2 = abs(abs(c) - 1) < ε || abs(λ0 - λ1) < ε ? (λ0 + λ1) / 2 : Math.atan2(b, a),\n\t          p = project(λ2, φ2),\n\t          x2 = p[0],\n\t          y2 = p[1],\n\t          dx2 = x2 - x0,\n\t          dy2 = y2 - y0,\n\t          dz = dy * dx2 - dx * dy2;\n\t      if (dz * dz / d2 > δ2 // perpendicular projected distance\n\t          || abs((dx * dx2 + dy * dy2) / d2 - .5) > .3 // midpoint close to an end\n\t          || a0 * a1 + b0 * b1 + c0 * c1 < cosMinDistance) { // angular distance\n\t        resampleLineTo(x0, y0, λ0, a0, b0, c0, x2, y2, λ2, a /= m, b /= m, c, depth, stream);\n\t        stream.point(x2, y2);\n\t        resampleLineTo(x2, y2, λ2, a, b, c, x1, y1, λ1, a1, b1, c1, depth, stream);\n\t      }\n\t    }\n\t  }\n\n\t  resample.precision = function(_) {\n\t    if (!arguments.length) return Math.sqrt(δ2);\n\t    maxDepth = (δ2 = _ * _) > 0 && 16;\n\t    return resample;\n\t  };\n\n\t  return resample;\n\t}\n\tvar d3f_arraySlice = [].slice,\n\t    d3f_array = function(list) { return d3f_arraySlice.call(list); }; // conversion for NodeLists\n\n\td3f.geo.transform = function(methods) {\n\t  return {\n\t    stream: function(stream) {\n\t      var transform = new d3f_geo_transform(stream);\n\t      for (var k in methods) transform[k] = methods[k];\n\t      return transform;\n\t    }\n\t  };\n\t};\n\n\tfunction d3f_geo_transform(stream) {\n\t  this.stream = stream;\n\t}\n\n\td3f_geo_transform.prototype = {\n\t  point: function(x, y) { this.stream.point(x, y); },\n\t  sphere: function() { this.stream.sphere(); },\n\t  lineStart: function() { this.stream.lineStart(); },\n\t  lineEnd: function() { this.stream.lineEnd(); },\n\t  polygonStart: function() { this.stream.polygonStart(); },\n\t  polygonEnd: function() { this.stream.polygonEnd(); }\n\t};\n\n\tfunction d3f_geo_transformPoint(stream, point) {\n\t  return {\n\t    point: point,\n\t    sphere: function() { stream.sphere(); },\n\t    lineStart: function() { stream.lineStart(); },\n\t    lineEnd: function() { stream.lineEnd(); },\n\t    polygonStart: function() { stream.polygonStart(); },\n\t    polygonEnd: function() { stream.polygonEnd(); },\n\t  };\n\t}\n\n\td3f.geo.projection = d3f_geo_projection;\n\td3f.geo.projectionMutator = d3f_geo_projectionMutator;\n\n\tfunction d3f_geo_projection(project) {\n\t  return d3f_geo_projectionMutator(function() { return project; })();\n\t}\n\n\tfunction d3f_geo_projectionMutator(projectAt) {\n\t  var project,\n\t      rotate,\n\t      projectRotate,\n\t      projectResample = d3f_geo_resample(function(x, y) { x = project(x, y); return [x[0] * k + δx, δy - x[1] * k]; }),\n\t      k = 150, // scale\n\t      x = 480, y = 250, // translate\n\t      λ = 0, φ = 0, // center\n\t      δλ = 0, δφ = 0, δγ = 0, // rotate\n\t      δx, δy, // center\n\t      preclip = d3f_geo_clipAntimeridian,\n\t      postclip = d3f_identity,\n\t      clipAngle = null,\n\t      clipExtent = null,\n\t      stream;\n\n\t  function projection(point) {\n\t    point = projectRotate(point[0] * d3f_radians, point[1] * d3f_radians);\n\t    return [point[0] * k + δx, δy - point[1] * k];\n\t  }\n\n\t  function invert(point) {\n\t    point = projectRotate.invert((point[0] - δx) / k, (δy - point[1]) / k);\n\t    return point && [point[0] * d3f_degrees, point[1] * d3f_degrees];\n\t  }\n\n\t  projection.stream = function(output) {\n\t    if (stream) stream.valid = false;\n\t    stream = d3f_geo_projectionRadians(preclip(rotate, projectResample(postclip(output))));\n\t    stream.valid = true; // allow caching by d3f.geo.path\n\t    return stream;\n\t  };\n\n\t  projection.clipAngle = function(_) {\n\t    if (!arguments.length) return clipAngle;\n\t    preclip = _ == null ? (clipAngle = _, d3f_geo_clipAntimeridian) : d3f_geo_clipCircle((clipAngle = +_) * d3f_radians);\n\t    return invalidate();\n\t  };\n\n\t  projection.clipExtent = function(_) {\n\t    if (!arguments.length) return clipExtent;\n\t    clipExtent = _;\n\t    postclip = _ ? d3f_geo_clipExtent(_[0][0], _[0][1], _[1][0], _[1][1]) : d3f_identity;\n\t    return invalidate();\n\t  };\n\n\t  projection.scale = function(_) {\n\t    if (!arguments.length) return k;\n\t    k = +_;\n\t    return reset();\n\t  };\n\n\t  projection.translate = function(_) {\n\t    if (!arguments.length) return [x, y];\n\t    x = +_[0];\n\t    y = +_[1];\n\t    return reset();\n\t  };\n\n\t  projection.center = function(_) {\n\t    if (!arguments.length) return [λ * d3f_degrees, φ * d3f_degrees];\n\t    λ = _[0] % 360 * d3f_radians;\n\t    φ = _[1] % 360 * d3f_radians;\n\t    return reset();\n\t  };\n\n\t  projection.rotate = function(_) {\n\t    if (!arguments.length) return [δλ * d3f_degrees, δφ * d3f_degrees, δγ * d3f_degrees];\n\t    δλ = _[0] % 360 * d3f_radians;\n\t    δφ = _[1] % 360 * d3f_radians;\n\t    δγ = _.length > 2 ? _[2] % 360 * d3f_radians : 0;\n\t    return reset();\n\t  };\n\n\t  d3f.rebind(projection, projectResample, \"precision\");\n\n\t  function reset() {\n\t    projectRotate = d3f_geo_compose(rotate = d3f_geo_rotation(δλ, δφ, δγ), project);\n\t    var center = project(λ, φ);\n\t    δx = x - center[0] * k;\n\t    δy = y + center[1] * k;\n\t    return invalidate();\n\t  }\n\n\t  function invalidate() {\n\t    if (stream) stream.valid = false, stream = null;\n\t    return projection;\n\t  }\n\n\t  return function() {\n\t    project = projectAt.apply(this, arguments);\n\t    projection.invert = project.invert && invert;\n\t    return reset();\n\t  };\n\t}\n\n\tfunction d3f_geo_projectionRadians(stream) {\n\t  return d3f_geo_transformPoint(stream, function(x, y) {\n\t    stream.point(x * d3f_radians, y * d3f_radians);\n\t  });\n\t}\n\n\tfunction d3f_geo_conic(projectAt) {\n\t  var φ0 = 0,\n\t      φ1 = π / 3,\n\t      m = d3f_geo_projectionMutator(projectAt),\n\t      p = m(φ0, φ1);\n\n\t  p.parallels = function(_) {\n\t    if (!arguments.length) return [φ0 / π * 180, φ1 / π * 180];\n\t    return m(φ0 = _[0] * π / 180, φ1 = _[1] * π / 180);\n\t  };\n\n\t  return p;\n\t}\n\n\tfunction d3f_geo_conicEqualArea(φ0, φ1) {\n\t  var sinφ0 = Math.sin(φ0),\n\t      n = (sinφ0 + Math.sin(φ1)) / 2,\n\t      C = 1 + sinφ0 * (2 * n - sinφ0),\n\t      ρ0 = Math.sqrt(C) / n;\n\n\t  function forward(λ, φ) {\n\t    var ρ = Math.sqrt(C - 2 * n * Math.sin(φ)) / n;\n\t    return [\n\t      ρ * Math.sin(λ *= n),\n\t      ρ0 - ρ * Math.cos(λ)\n\t    ];\n\t  }\n\n\t  forward.invert = function(x, y) {\n\t    var ρ0_y = ρ0 - y;\n\t    return [\n\t      Math.atan2(x, ρ0_y) / n,\n\t      d3f_asin((C - (x * x + ρ0_y * ρ0_y) * n * n) / (2 * n))\n\t    ];\n\t  };\n\n\t  return forward;\n\t}\n\n\t(d3f.geo.conicEqualArea = function() {\n\t  return d3f_geo_conic(d3f_geo_conicEqualArea);\n\t}).raw = d3f_geo_conicEqualArea;\n\n\t// ESRI:102003\n\td3f.geo.albers = function() {\n\t  return d3f.geo.conicEqualArea()\n\t      .rotate([96, 0])\n\t      .center([-.6, 38.7])\n\t      .parallels([29.5, 45.5])\n\t      .scale(1070);\n\t};\n\n\t// A composite projection for the United States, configured by default for\n\t// 960×500. Also works quite well at 960×600 with scale 1285. The set of\n\t// standard parallels for each region comes from USGS, which is published here:\n\t// http://egsc.usgs.gov/isb/pubs/MapProjections/projections.html#albers\n\td3f.geo.albersUsa = function() {\n\t  var lower48 = d3f.geo.albers();\n\n\t  // EPSG:3338\n\t  var alaska = d3f.geo.conicEqualArea()\n\t      .rotate([154, 0])\n\t      .center([-2, 58.5])\n\t      .parallels([55, 65]);\n\n\t  // ESRI:102007\n\t  var hawaii = d3f.geo.conicEqualArea()\n\t      .rotate([157, 0])\n\t      .center([-3, 19.9])\n\t      .parallels([8, 18]);\n\n\t  var point,\n\t      pointStream = {point: function(x, y) { point = [x, y]; }},\n\t      lower48Point,\n\t      alaskaPoint,\n\t      hawaiiPoint;\n\n\t  function albersUsa(coordinates) {\n\t    var x = coordinates[0], y = coordinates[1];\n\t    point = null;\n\t    (lower48Point(x, y), point)\n\t        || (alaskaPoint(x, y), point)\n\t        || hawaiiPoint(x, y);\n\t    return point;\n\t  }\n\n\t  albersUsa.invert = function(coordinates) {\n\t    var k = lower48.scale(),\n\t        t = lower48.translate(),\n\t        x = (coordinates[0] - t[0]) / k,\n\t        y = (coordinates[1] - t[1]) / k;\n\t    return (y >= .120 && y < .234 && x >= -.425 && x < -.214 ? alaska\n\t        : y >= .166 && y < .234 && x >= -.214 && x < -.115 ? hawaii\n\t        : lower48).invert(coordinates);\n\t  };\n\n\t  // A naïve multi-projection stream.\n\t  // The projections must have mutually exclusive clip regions on the sphere,\n\t  // as this will avoid emitting interleaving lines and polygons.\n\t  albersUsa.stream = function(stream) {\n\t    var lower48Stream = lower48.stream(stream),\n\t        alaskaStream = alaska.stream(stream),\n\t        hawaiiStream = hawaii.stream(stream);\n\t    return {\n\t      point: function(x, y) {\n\t        lower48Stream.point(x, y);\n\t        alaskaStream.point(x, y);\n\t        hawaiiStream.point(x, y);\n\t      },\n\t      sphere: function() {\n\t        lower48Stream.sphere();\n\t        alaskaStream.sphere();\n\t        hawaiiStream.sphere();\n\t      },\n\t      lineStart: function() {\n\t        lower48Stream.lineStart();\n\t        alaskaStream.lineStart();\n\t        hawaiiStream.lineStart();\n\t      },\n\t      lineEnd: function() {\n\t        lower48Stream.lineEnd();\n\t        alaskaStream.lineEnd();\n\t        hawaiiStream.lineEnd();\n\t      },\n\t      polygonStart: function() {\n\t        lower48Stream.polygonStart();\n\t        alaskaStream.polygonStart();\n\t        hawaiiStream.polygonStart();\n\t      },\n\t      polygonEnd: function() {\n\t        lower48Stream.polygonEnd();\n\t        alaskaStream.polygonEnd();\n\t        hawaiiStream.polygonEnd();\n\t      }\n\t    };\n\t  };\n\n\t  albersUsa.precision = function(_) {\n\t    if (!arguments.length) return lower48.precision();\n\t    lower48.precision(_);\n\t    alaska.precision(_);\n\t    hawaii.precision(_);\n\t    return albersUsa;\n\t  };\n\n\t  albersUsa.scale = function(_) {\n\t    if (!arguments.length) return lower48.scale();\n\t    lower48.scale(_);\n\t    alaska.scale(_ * .35);\n\t    hawaii.scale(_);\n\t    return albersUsa.translate(lower48.translate());\n\t  };\n\n\t  albersUsa.translate = function(_) {\n\t    if (!arguments.length) return lower48.translate();\n\t    var k = lower48.scale(), x = +_[0], y = +_[1];\n\n\t    lower48Point = lower48\n\t        .translate(_)\n\t        .clipExtent([[x - .455 * k, y - .238 * k], [x + .455 * k, y + .238 * k]])\n\t        .stream(pointStream).point;\n\n\t    alaskaPoint = alaska\n\t        .translate([x - .307 * k, y + .201 * k])\n\t        .clipExtent([[x - .425 * k + ε, y + .120 * k + ε], [x - .214 * k - ε, y + .234 * k - ε]])\n\t        .stream(pointStream).point;\n\n\t    hawaiiPoint = hawaii\n\t        .translate([x - .205 * k, y + .212 * k])\n\t        .clipExtent([[x - .214 * k + ε, y + .166 * k + ε], [x - .115 * k - ε, y + .234 * k - ε]])\n\t        .stream(pointStream).point;\n\n\t    return albersUsa;\n\t  };\n\n\t  return albersUsa.scale(1070);\n\t};\n\n\td3f.geo.bounds = (function() {\n\t  var λ0, φ0, λ1, φ1, // bounds\n\t      λ_, // previous λ-coordinate\n\t      λ__, φ__, // first point\n\t      p0, // previous 3D point\n\t      dλSum,\n\t      ranges,\n\t      range;\n\n\t  var bound = {\n\t    point: point,\n\t    lineStart: lineStart,\n\t    lineEnd: lineEnd,\n\n\t    polygonStart: function() {\n\t      bound.point = ringPoint;\n\t      bound.lineStart = ringStart;\n\t      bound.lineEnd = ringEnd;\n\t      dλSum = 0;\n\t      d3f_geo_area.polygonStart();\n\t    },\n\t    polygonEnd: function() {\n\t      d3f_geo_area.polygonEnd();\n\t      bound.point = point;\n\t      bound.lineStart = lineStart;\n\t      bound.lineEnd = lineEnd;\n\t      if (d3f_geo_areaRingSum < 0) λ0 = -(λ1 = 180), φ0 = -(φ1 = 90);\n\t      else if (dλSum > ε) φ1 = 90;\n\t      else if (dλSum < -ε) φ0 = -90;\n\t      range[0] = λ0, range[1] = λ1;\n\t    }\n\t  };\n\n\t  function point(λ, φ) {\n\t    ranges.push(range = [λ0 = λ, λ1 = λ]);\n\t    if (φ < φ0) φ0 = φ;\n\t    if (φ > φ1) φ1 = φ;\n\t  }\n\n\t  function linePoint(λ, φ) {\n\t    var p = d3f_geo_cartesian([λ * d3f_radians, φ * d3f_radians]);\n\t    if (p0) {\n\t      var normal = d3f_geo_cartesianCross(p0, p),\n\t          equatorial = [normal[1], -normal[0], 0],\n\t          inflection = d3f_geo_cartesianCross(equatorial, normal);\n\t      d3f_geo_cartesianNormalize(inflection);\n\t      inflection = d3f_geo_spherical(inflection);\n\t      var dλ = λ - λ_,\n\t          s = dλ > 0 ? 1 : -1,\n\t          λi = inflection[0] * d3f_degrees * s,\n\t          antimeridian = abs(dλ) > 180;\n\t      if (antimeridian ^ (s * λ_ < λi && λi < s * λ)) {\n\t        var φi = inflection[1] * d3f_degrees;\n\t        if (φi > φ1) φ1 = φi;\n\t      } else if (λi = (λi + 360) % 360 - 180, antimeridian ^ (s * λ_ < λi && λi < s * λ)) {\n\t        var φi = -inflection[1] * d3f_degrees;\n\t        if (φi < φ0) φ0 = φi;\n\t      } else {\n\t        if (φ < φ0) φ0 = φ;\n\t        if (φ > φ1) φ1 = φ;\n\t      }\n\t      if (antimeridian) {\n\t        if (λ < λ_) {\n\t          if (angle(λ0, λ) > angle(λ0, λ1)) λ1 = λ;\n\t        } else {\n\t          if (angle(λ, λ1) > angle(λ0, λ1)) λ0 = λ;\n\t        }\n\t      } else {\n\t        if (λ1 >= λ0) {\n\t          if (λ < λ0) λ0 = λ;\n\t          if (λ > λ1) λ1 = λ;\n\t        } else {\n\t          if (λ > λ_) {\n\t            if (angle(λ0, λ) > angle(λ0, λ1)) λ1 = λ;\n\t          } else {\n\t            if (angle(λ, λ1) > angle(λ0, λ1)) λ0 = λ;\n\t          }\n\t        }\n\t      }\n\t    } else {\n\t      point(λ, φ);\n\t    }\n\t    p0 = p, λ_ = λ;\n\t  }\n\n\t  function lineStart() { bound.point = linePoint; }\n\t  function lineEnd() {\n\t    range[0] = λ0, range[1] = λ1;\n\t    bound.point = point;\n\t    p0 = null;\n\t  }\n\n\t  function ringPoint(λ, φ) {\n\t    if (p0) {\n\t      var dλ = λ - λ_;\n\t      dλSum += abs(dλ) > 180 ? dλ + (dλ > 0 ? 360 : -360) : dλ;\n\t    } else λ__ = λ, φ__ = φ;\n\t    d3f_geo_area.point(λ, φ);\n\t    linePoint(λ, φ);\n\t  }\n\n\t  function ringStart() {\n\t    d3f_geo_area.lineStart();\n\t  }\n\n\t  function ringEnd() {\n\t    ringPoint(λ__, φ__);\n\t    d3f_geo_area.lineEnd();\n\t    if (abs(dλSum) > ε) λ0 = -(λ1 = 180);\n\t    range[0] = λ0, range[1] = λ1;\n\t    p0 = null;\n\t  }\n\n\t  // Finds the left-right distance between two longitudes.\n\t  // This is almost the same as (λ1 - λ0 + 360°) % 360°, except that we want\n\t  // the distance between ±180° to be 360°.\n\t  function angle(λ0, λ1) { return (λ1 -= λ0) < 0 ? λ1 + 360 : λ1; }\n\n\t  function compareRanges(a, b) { return a[0] - b[0]; }\n\n\t  function withinRange(x, range) {\n\t    return range[0] <= range[1] ? range[0] <= x && x <= range[1] : x < range[0] || range[1] < x;\n\t  }\n\n\t  return function(feature) {\n\t    φ1 = λ1 = -(λ0 = φ0 = Infinity);\n\t    ranges = [];\n\n\t    d3f.geo.stream(feature, bound);\n\n\t    var n = ranges.length;\n\t    if (n) {\n\t      // First, sort ranges by their minimum longitudes.\n\t      ranges.sort(compareRanges);\n\n\t      // Then, merge any ranges that overlap.\n\t      for (var i = 1, a = ranges[0], b, merged = [a]; i < n; ++i) {\n\t        b = ranges[i];\n\t        if (withinRange(b[0], a) || withinRange(b[1], a)) {\n\t          if (angle(a[0], b[1]) > angle(a[0], a[1])) a[1] = b[1];\n\t          if (angle(b[0], a[1]) > angle(a[0], a[1])) a[0] = b[0];\n\t        } else {\n\t          merged.push(a = b);\n\t        }\n\t      }\n\n\t      // Finally, find the largest gap between the merged ranges.\n\t      // The final bounding box will be the inverse of this gap.\n\t      var best = -Infinity, dλ;\n\t      for (var n = merged.length - 1, i = 0, a = merged[n], b; i <= n; a = b, ++i) {\n\t        b = merged[i];\n\t        if ((dλ = angle(a[1], b[0])) > best) best = dλ, λ0 = b[0], λ1 = a[1];\n\t      }\n\t    }\n\t    ranges = range = null;\n\n\t    return λ0 === Infinity || φ0 === Infinity\n\t        ? [[NaN, NaN], [NaN, NaN]]\n\t        : [[λ0, φ0], [λ1, φ1]];\n\t  };\n\t})();\n\n\td3f.geo.centroid = function(object) {\n\t  d3f_geo_centroidW0 = d3f_geo_centroidW1 =\n\t  d3f_geo_centroidX0 = d3f_geo_centroidY0 = d3f_geo_centroidZ0 =\n\t  d3f_geo_centroidX1 = d3f_geo_centroidY1 = d3f_geo_centroidZ1 =\n\t  d3f_geo_centroidX2 = d3f_geo_centroidY2 = d3f_geo_centroidZ2 = 0;\n\t  d3f.geo.stream(object, d3f_geo_centroid);\n\n\t  var x = d3f_geo_centroidX2,\n\t      y = d3f_geo_centroidY2,\n\t      z = d3f_geo_centroidZ2,\n\t      m = x * x + y * y + z * z;\n\n\t  // If the area-weighted centroid is undefined, fall back to length-weighted centroid.\n\t  if (m < ε2) {\n\t    x = d3f_geo_centroidX1, y = d3f_geo_centroidY1, z = d3f_geo_centroidZ1;\n\t    // If the feature has zero length, fall back to arithmetic mean of point vectors.\n\t    if (d3f_geo_centroidW1 < ε) x = d3f_geo_centroidX0, y = d3f_geo_centroidY0, z = d3f_geo_centroidZ0;\n\t    m = x * x + y * y + z * z;\n\t    // If the feature still has an undefined centroid, then return.\n\t    if (m < ε2) return [NaN, NaN];\n\t  }\n\n\t  return [Math.atan2(y, x) * d3f_degrees, d3f_asin(z / Math.sqrt(m)) * d3f_degrees];\n\t};\n\n\tvar d3f_geo_centroidW0,\n\t    d3f_geo_centroidW1,\n\t    d3f_geo_centroidX0,\n\t    d3f_geo_centroidY0,\n\t    d3f_geo_centroidZ0,\n\t    d3f_geo_centroidX1,\n\t    d3f_geo_centroidY1,\n\t    d3f_geo_centroidZ1,\n\t    d3f_geo_centroidX2,\n\t    d3f_geo_centroidY2,\n\t    d3f_geo_centroidZ2;\n\n\tvar d3f_geo_centroid = {\n\t  sphere: d3f_noop,\n\t  point: d3f_geo_centroidPoint,\n\t  lineStart: d3f_geo_centroidLineStart,\n\t  lineEnd: d3f_geo_centroidLineEnd,\n\t  polygonStart: function() {\n\t    d3f_geo_centroid.lineStart = d3f_geo_centroidRingStart;\n\t  },\n\t  polygonEnd: function() {\n\t    d3f_geo_centroid.lineStart = d3f_geo_centroidLineStart;\n\t  }\n\t};\n\n\t// Arithmetic mean of Cartesian vectors.\n\tfunction d3f_geo_centroidPoint(λ, φ) {\n\t  λ *= d3f_radians;\n\t  var cosφ = Math.cos(φ *= d3f_radians);\n\t  d3f_geo_centroidPointXYZ(cosφ * Math.cos(λ), cosφ * Math.sin(λ), Math.sin(φ));\n\t}\n\n\tfunction d3f_geo_centroidPointXYZ(x, y, z) {\n\t  ++d3f_geo_centroidW0;\n\t  d3f_geo_centroidX0 += (x - d3f_geo_centroidX0) / d3f_geo_centroidW0;\n\t  d3f_geo_centroidY0 += (y - d3f_geo_centroidY0) / d3f_geo_centroidW0;\n\t  d3f_geo_centroidZ0 += (z - d3f_geo_centroidZ0) / d3f_geo_centroidW0;\n\t}\n\n\tfunction d3f_geo_centroidLineStart() {\n\t  var x0, y0, z0; // previous point\n\n\t  d3f_geo_centroid.point = function(λ, φ) {\n\t    λ *= d3f_radians;\n\t    var cosφ = Math.cos(φ *= d3f_radians);\n\t    x0 = cosφ * Math.cos(λ);\n\t    y0 = cosφ * Math.sin(λ);\n\t    z0 = Math.sin(φ);\n\t    d3f_geo_centroid.point = nextPoint;\n\t    d3f_geo_centroidPointXYZ(x0, y0, z0);\n\t  };\n\n\t  function nextPoint(λ, φ) {\n\t    λ *= d3f_radians;\n\t    var cosφ = Math.cos(φ *= d3f_radians),\n\t        x = cosφ * Math.cos(λ),\n\t        y = cosφ * Math.sin(λ),\n\t        z = Math.sin(φ),\n\t        w = Math.atan2(\n\t          Math.sqrt((w = y0 * z - z0 * y) * w + (w = z0 * x - x0 * z) * w + (w = x0 * y - y0 * x) * w),\n\t          x0 * x + y0 * y + z0 * z);\n\t    d3f_geo_centroidW1 += w;\n\t    d3f_geo_centroidX1 += w * (x0 + (x0 = x));\n\t    d3f_geo_centroidY1 += w * (y0 + (y0 = y));\n\t    d3f_geo_centroidZ1 += w * (z0 + (z0 = z));\n\t    d3f_geo_centroidPointXYZ(x0, y0, z0);\n\t  }\n\t}\n\n\tfunction d3f_geo_centroidLineEnd() {\n\t  d3f_geo_centroid.point = d3f_geo_centroidPoint;\n\t}\n\n\t// See J. E. Brock, The Inertia Tensor for a Spherical Triangle,\n\t// J. Applied Mechanics 42, 239 (1975).\n\tfunction d3f_geo_centroidRingStart() {\n\t  var λ00, φ00, // first point\n\t      x0, y0, z0; // previous point\n\n\t  d3f_geo_centroid.point = function(λ, φ) {\n\t    λ00 = λ, φ00 = φ;\n\t    d3f_geo_centroid.point = nextPoint;\n\t    λ *= d3f_radians;\n\t    var cosφ = Math.cos(φ *= d3f_radians);\n\t    x0 = cosφ * Math.cos(λ);\n\t    y0 = cosφ * Math.sin(λ);\n\t    z0 = Math.sin(φ);\n\t    d3f_geo_centroidPointXYZ(x0, y0, z0);\n\t  };\n\n\t  d3f_geo_centroid.lineEnd = function() {\n\t    nextPoint(λ00, φ00);\n\t    d3f_geo_centroid.lineEnd = d3f_geo_centroidLineEnd;\n\t    d3f_geo_centroid.point = d3f_geo_centroidPoint;\n\t  };\n\n\t  function nextPoint(λ, φ) {\n\t    λ *= d3f_radians;\n\t    var cosφ = Math.cos(φ *= d3f_radians),\n\t        x = cosφ * Math.cos(λ),\n\t        y = cosφ * Math.sin(λ),\n\t        z = Math.sin(φ),\n\t        cx = y0 * z - z0 * y,\n\t        cy = z0 * x - x0 * z,\n\t        cz = x0 * y - y0 * x,\n\t        m = Math.sqrt(cx * cx + cy * cy + cz * cz),\n\t        u = x0 * x + y0 * y + z0 * z,\n\t        v = m && -d3f_acos(u) / m, // area weight\n\t        w = Math.atan2(m, u); // line weight\n\t    d3f_geo_centroidX2 += v * cx;\n\t    d3f_geo_centroidY2 += v * cy;\n\t    d3f_geo_centroidZ2 += v * cz;\n\t    d3f_geo_centroidW1 += w;\n\t    d3f_geo_centroidX1 += w * (x0 + (x0 = x));\n\t    d3f_geo_centroidY1 += w * (y0 + (y0 = y));\n\t    d3f_geo_centroidZ1 += w * (z0 + (z0 = z));\n\t    d3f_geo_centroidPointXYZ(x0, y0, z0);\n\t  }\n\t}\n\n\t// TODO Unify this code with d3f.geom.polygon area?\n\n\tvar d3f_geo_pathAreaSum, d3f_geo_pathAreaPolygon, d3f_geo_pathArea = {\n\t  point: d3f_noop,\n\t  lineStart: d3f_noop,\n\t  lineEnd: d3f_noop,\n\n\t  // Only count area for polygon rings.\n\t  polygonStart: function() {\n\t    d3f_geo_pathAreaPolygon = 0;\n\t    d3f_geo_pathArea.lineStart = d3f_geo_pathAreaRingStart;\n\t  },\n\t  polygonEnd: function() {\n\t    d3f_geo_pathArea.lineStart = d3f_geo_pathArea.lineEnd = d3f_geo_pathArea.point = d3f_noop;\n\t    d3f_geo_pathAreaSum += abs(d3f_geo_pathAreaPolygon / 2);\n\t  }\n\t};\n\n\tfunction d3f_geo_pathAreaRingStart() {\n\t  var x00, y00, x0, y0;\n\n\t  // For the first point, …\n\t  d3f_geo_pathArea.point = function(x, y) {\n\t    d3f_geo_pathArea.point = nextPoint;\n\t    x00 = x0 = x, y00 = y0 = y;\n\t  };\n\n\t  // For subsequent points, …\n\t  function nextPoint(x, y) {\n\t    d3f_geo_pathAreaPolygon += y0 * x - x0 * y;\n\t    x0 = x, y0 = y;\n\t  }\n\n\t  // For the last point, return to the start.\n\t  d3f_geo_pathArea.lineEnd = function() {\n\t    nextPoint(x00, y00);\n\t  };\n\t}\n\n\tvar d3f_geo_pathBoundsX0,\n\t    d3f_geo_pathBoundsY0,\n\t    d3f_geo_pathBoundsX1,\n\t    d3f_geo_pathBoundsY1;\n\n\tvar d3f_geo_pathBounds = {\n\t  point: d3f_geo_pathBoundsPoint,\n\t  lineStart: d3f_noop,\n\t  lineEnd: d3f_noop,\n\t  polygonStart: d3f_noop,\n\t  polygonEnd: d3f_noop\n\t};\n\n\tfunction d3f_geo_pathBoundsPoint(x, y) {\n\t  if (x < d3f_geo_pathBoundsX0) d3f_geo_pathBoundsX0 = x;\n\t  if (x > d3f_geo_pathBoundsX1) d3f_geo_pathBoundsX1 = x;\n\t  if (y < d3f_geo_pathBoundsY0) d3f_geo_pathBoundsY0 = y;\n\t  if (y > d3f_geo_pathBoundsY1) d3f_geo_pathBoundsY1 = y;\n\t}\n\tfunction d3f_geo_pathBuffer() {\n\t  var pointCircle = d3f_geo_pathBufferCircle(4.5),\n\t      buffer = [];\n\n\t  var stream = {\n\t    point: point,\n\n\t    // While inside a line, override point to moveTo then lineTo.\n\t    lineStart: function() { stream.point = pointLineStart; },\n\t    lineEnd: lineEnd,\n\n\t    // While inside a polygon, override lineEnd to closePath.\n\t    polygonStart: function() { stream.lineEnd = lineEndPolygon; },\n\t    polygonEnd: function() { stream.lineEnd = lineEnd; stream.point = point; },\n\n\t    pointRadius: function(_) {\n\t      pointCircle = d3f_geo_pathBufferCircle(_);\n\t      return stream;\n\t    },\n\n\t    result: function() {\n\t      if (buffer.length) {\n\t        var result = buffer.join(\"\");\n\t        buffer = [];\n\t        return result;\n\t      }\n\t    }\n\t  };\n\n\t  function point(x, y) {\n\t    buffer.push(\"M\", x, \",\", y, pointCircle);\n\t  }\n\n\t  function pointLineStart(x, y) {\n\t    buffer.push(\"M\", x, \",\", y);\n\t    stream.point = pointLine;\n\t  }\n\n\t  function pointLine(x, y) {\n\t    buffer.push(\"L\", x, \",\", y);\n\t  }\n\n\t  function lineEnd() {\n\t    stream.point = point;\n\t  }\n\n\t  function lineEndPolygon() {\n\t    buffer.push(\"Z\");\n\t  }\n\n\t  return stream;\n\t}\n\n\tfunction d3f_geo_pathBufferCircle(radius) {\n\t  return \"m0,\" + radius\n\t      + \"a\" + radius + \",\" + radius + \" 0 1,1 0,\" + -2 * radius\n\t      + \"a\" + radius + \",\" + radius + \" 0 1,1 0,\" + 2 * radius\n\t      + \"z\";\n\t}\n\n\t// TODO Unify this code with d3f.geom.polygon centroid?\n\t// TODO Enforce positive area for exterior, negative area for interior?\n\n\tvar d3f_geo_pathCentroid = {\n\t  point: d3f_geo_pathCentroidPoint,\n\n\t  // For lines, weight by length.\n\t  lineStart: d3f_geo_pathCentroidLineStart,\n\t  lineEnd: d3f_geo_pathCentroidLineEnd,\n\n\t  // For polygons, weight by area.\n\t  polygonStart: function() {\n\t    d3f_geo_pathCentroid.lineStart = d3f_geo_pathCentroidRingStart;\n\t  },\n\t  polygonEnd: function() {\n\t    d3f_geo_pathCentroid.point = d3f_geo_pathCentroidPoint;\n\t    d3f_geo_pathCentroid.lineStart = d3f_geo_pathCentroidLineStart;\n\t    d3f_geo_pathCentroid.lineEnd = d3f_geo_pathCentroidLineEnd;\n\t  }\n\t};\n\n\tfunction d3f_geo_pathCentroidPoint(x, y) {\n\t  d3f_geo_centroidX0 += x;\n\t  d3f_geo_centroidY0 += y;\n\t  ++d3f_geo_centroidZ0;\n\t}\n\n\tfunction d3f_geo_pathCentroidLineStart() {\n\t  var x0, y0;\n\n\t  d3f_geo_pathCentroid.point = function(x, y) {\n\t    d3f_geo_pathCentroid.point = nextPoint;\n\t    d3f_geo_pathCentroidPoint(x0 = x, y0 = y);\n\t  };\n\n\t  function nextPoint(x, y) {\n\t    var dx = x - x0, dy = y - y0, z = Math.sqrt(dx * dx + dy * dy);\n\t    d3f_geo_centroidX1 += z * (x0 + x) / 2;\n\t    d3f_geo_centroidY1 += z * (y0 + y) / 2;\n\t    d3f_geo_centroidZ1 += z;\n\t    d3f_geo_pathCentroidPoint(x0 = x, y0 = y);\n\t  }\n\t}\n\n\tfunction d3f_geo_pathCentroidLineEnd() {\n\t  d3f_geo_pathCentroid.point = d3f_geo_pathCentroidPoint;\n\t}\n\n\tfunction d3f_geo_pathCentroidRingStart() {\n\t  var x00, y00, x0, y0;\n\n\t  // For the first point, …\n\t  d3f_geo_pathCentroid.point = function(x, y) {\n\t    d3f_geo_pathCentroid.point = nextPoint;\n\t    d3f_geo_pathCentroidPoint(x00 = x0 = x, y00 = y0 = y);\n\t  };\n\n\t  // For subsequent points, …\n\t  function nextPoint(x, y) {\n\t    var dx = x - x0, dy = y - y0, z = Math.sqrt(dx * dx + dy * dy);\n\t    d3f_geo_centroidX1 += z * (x0 + x) / 2;\n\t    d3f_geo_centroidY1 += z * (y0 + y) / 2;\n\t    d3f_geo_centroidZ1 += z;\n\n\t    z = y0 * x - x0 * y;\n\t    d3f_geo_centroidX2 += z * (x0 + x);\n\t    d3f_geo_centroidY2 += z * (y0 + y);\n\t    d3f_geo_centroidZ2 += z * 3;\n\t    d3f_geo_pathCentroidPoint(x0 = x, y0 = y);\n\t  }\n\n\t  // For the last point, return to the start.\n\t  d3f_geo_pathCentroid.lineEnd = function() {\n\t    nextPoint(x00, y00);\n\t  };\n\t}\n\n\tfunction d3f_geo_pathContext(context) {\n\t  var pointRadius = 4.5;\n\n\t  var stream = {\n\t    point: point,\n\n\t    // While inside a line, override point to moveTo then lineTo.\n\t    lineStart: function() { stream.point = pointLineStart; },\n\t    lineEnd: lineEnd,\n\n\t    // While inside a polygon, override lineEnd to closePath.\n\t    polygonStart: function() { stream.lineEnd = lineEndPolygon; },\n\t    polygonEnd: function() { stream.lineEnd = lineEnd; stream.point = point; },\n\n\t    pointRadius: function(_) {\n\t      pointRadius = _;\n\t      return stream;\n\t    },\n\n\t    result: d3f_noop\n\t  };\n\n\t  function point(x, y) {\n\t    context.moveTo(x + pointRadius, y);\n\t    context.arc(x, y, pointRadius, 0, τ);\n\t  }\n\n\t  function pointLineStart(x, y) {\n\t    context.moveTo(x, y);\n\t    stream.point = pointLine;\n\t  }\n\n\t  function pointLine(x, y) {\n\t    context.lineTo(x, y);\n\t  }\n\n\t  function lineEnd() {\n\t    stream.point = point;\n\t  }\n\n\t  function lineEndPolygon() {\n\t    context.closePath();\n\t  }\n\n\t  return stream;\n\t}\n\n\td3f.geo.path = function() {\n\t  var pointRadius = 4.5,\n\t      projection,\n\t      context,\n\t      projectStream,\n\t      contextStream,\n\t      cacheStream;\n\n\t  function path(object) {\n\t    if (object) {\n\t      if (typeof pointRadius === \"function\") contextStream.pointRadius(+pointRadius.apply(this, arguments));\n\t      if (!cacheStream || !cacheStream.valid) cacheStream = projectStream(contextStream);\n\t      d3f.geo.stream(object, cacheStream);\n\t    }\n\t    return contextStream.result();\n\t  }\n\n\t  path.area = function(object) {\n\t    d3f_geo_pathAreaSum = 0;\n\t    d3f.geo.stream(object, projectStream(d3f_geo_pathArea));\n\t    return d3f_geo_pathAreaSum;\n\t  };\n\n\t  path.centroid = function(object) {\n\t    d3f_geo_centroidX0 = d3f_geo_centroidY0 = d3f_geo_centroidZ0 =\n\t    d3f_geo_centroidX1 = d3f_geo_centroidY1 = d3f_geo_centroidZ1 =\n\t    d3f_geo_centroidX2 = d3f_geo_centroidY2 = d3f_geo_centroidZ2 = 0;\n\t    d3f.geo.stream(object, projectStream(d3f_geo_pathCentroid));\n\t    return d3f_geo_centroidZ2 ? [d3f_geo_centroidX2 / d3f_geo_centroidZ2, d3f_geo_centroidY2 / d3f_geo_centroidZ2]\n\t        : d3f_geo_centroidZ1 ? [d3f_geo_centroidX1 / d3f_geo_centroidZ1, d3f_geo_centroidY1 / d3f_geo_centroidZ1]\n\t        : d3f_geo_centroidZ0 ? [d3f_geo_centroidX0 / d3f_geo_centroidZ0, d3f_geo_centroidY0 / d3f_geo_centroidZ0]\n\t        : [NaN, NaN];\n\t  };\n\n\t  path.bounds = function(object) {\n\t    d3f_geo_pathBoundsX1 = d3f_geo_pathBoundsY1 = -(d3f_geo_pathBoundsX0 = d3f_geo_pathBoundsY0 = Infinity);\n\t    d3f.geo.stream(object, projectStream(d3f_geo_pathBounds));\n\t    return [[d3f_geo_pathBoundsX0, d3f_geo_pathBoundsY0], [d3f_geo_pathBoundsX1, d3f_geo_pathBoundsY1]];\n\t  };\n\n\t  path.projection = function(_) {\n\t    if (!arguments.length) return projection;\n\t    projectStream = (projection = _) ? _.stream || d3f_geo_pathProjectStream(_) : d3f_identity;\n\t    return reset();\n\t  };\n\n\t  path.context = function(_) {\n\t    if (!arguments.length) return context;\n\t    contextStream = (context = _) == null ? new d3f_geo_pathBuffer : new d3f_geo_pathContext(_);\n\t    if (typeof pointRadius !== \"function\") contextStream.pointRadius(pointRadius);\n\t    return reset();\n\t  };\n\n\t  path.pointRadius = function(_) {\n\t    if (!arguments.length) return pointRadius;\n\t    pointRadius = typeof _ === \"function\" ? _ : (contextStream.pointRadius(+_), +_);\n\t    return path;\n\t  };\n\n\t  function reset() {\n\t    cacheStream = null;\n\t    return path;\n\t  }\n\n\t  return path.projection(d3f.geo.albersUsa()).context(null);\n\t};\n\n\tfunction d3f_geo_pathProjectStream(project) {\n\t  var resample = d3f_geo_resample(function(x, y) { return project([x * d3f_degrees, y * d3f_degrees]); });\n\t  return function(stream) { return d3f_geo_projectionRadians(resample(stream)); };\n\t}\n\n\tfunction d3f_geo_mercator(λ, φ) {\n\t  return [λ, Math.log(Math.tan(π / 4 + φ / 2))];\n\t}\n\n\td3f_geo_mercator.invert = function(x, y) {\n\t  return [x, 2 * Math.atan(Math.exp(y)) - halfπ];\n\t};\n\n\tfunction d3f_geo_mercatorProjection(project) {\n\t  var m = d3f_geo_projection(project),\n\t      scale = m.scale,\n\t      translate = m.translate,\n\t      clipExtent = m.clipExtent,\n\t      clipAuto;\n\n\t  m.scale = function() {\n\t    var v = scale.apply(m, arguments);\n\t    return v === m ? (clipAuto ? m.clipExtent(null) : m) : v;\n\t  };\n\n\t  m.translate = function() {\n\t    var v = translate.apply(m, arguments);\n\t    return v === m ? (clipAuto ? m.clipExtent(null) : m) : v;\n\t  };\n\n\t  m.clipExtent = function(_) {\n\t    var v = clipExtent.apply(m, arguments);\n\t    if (v === m) {\n\t      if (clipAuto = _ == null) {\n\t        var k = π * scale(), t = translate();\n\t        clipExtent([[t[0] - k, t[1] - k], [t[0] + k, t[1] + k]]);\n\t      }\n\t    } else if (clipAuto) {\n\t      v = null;\n\t    }\n\t    return v;\n\t  };\n\n\t  return m.clipExtent(null);\n\t}\n\n\t(d3f.geo.mercator = function() {\n\t  return d3f_geo_mercatorProjection(d3f_geo_mercator);\n\t}).raw = d3f_geo_mercator;\n\t  if (true) !(__WEBPACK_AMD_DEFINE_FACTORY__ = (d3f), __WEBPACK_AMD_DEFINE_RESULT__ = (typeof __WEBPACK_AMD_DEFINE_FACTORY__ === 'function' ? (__WEBPACK_AMD_DEFINE_FACTORY__.call(exports, __webpack_require__, exports, module)) : __WEBPACK_AMD_DEFINE_FACTORY__), __WEBPACK_AMD_DEFINE_RESULT__ !== undefined && (module.exports = __WEBPACK_AMD_DEFINE_RESULT__));\n\t  else if (typeof module === \"object\" && module.exports) module.exports = d3f;\n\t  this.d3f = d3f;\n\t}();\n\n\n/***/ },\n/* 2 */\n/***/ function(module, exports) {\n\n\tvar Helpers = {}\n\n\tHelpers.copy = function(o) {\n\t  return (o instanceof Array)\n\t  ? o.map(Helpers.copy)\n\t  : (typeof o === \"string\" || typeof o === \"number\")\n\t  ? o\n\t  : copyObject(o);\n\t}\n\n\tfunction copyObject(o) {\n\t  var obj = {};\n\t  for (var k in o) obj[k] = Helpers.copy(o[k]);\n\t  return obj;\n\t}\n\n\t//Grouping cosArctan and sinArctan (see below)\n\tHelpers.arctans = function(dx, dy){\n\t  var div = dx/dy,\n\t      sqrt = Math.sqrt(1+(div*div)),\n\t      signedSqrt = (dy > 0) ? sqrt : -sqrt,\n\t      cos = 1 / signedSqrt,\n\t      sin = div * cos;\n\n\t  return {\n\t    cos: cos,\n\t    sin: sin\n\t  }\n\t}\n\n\t/*\n\n\tfunction cosArctan(dx,dy){\n\t  var div = dx/dy;\n\t  return (dy>0)?\n\t  (1/Math.sqrt(1+(div*div))):\n\t  (-1/Math.sqrt(1+(div*div)));\n\t}\n\n\tfunction sinArctan(dx,dy){\n\t  var div = dx/dy;\n\t  return (dy>0)?\n\t  (div/Math.sqrt(1+(div*div))):\n\t  (-div/Math.sqrt(1+(div*div)));\n\t}\n\n\t*/\n\n\t// TODO what's this method for ? To optimize\n\tHelpers.object = function(arcs, o) {\n\t  function arc(i, points) {\n\t    if (points.length) points.pop();\n\t    for (var a = arcs[i < 0 ? ~i : i], k = 0, n = a.length; k < n; ++k) {\n\t      points.push(a[k]);\n\t    }\n\t    if (i < 0) reverse(points, n);\n\t  }\n\n\t  function line(arcs) {\n\t    var points = [];\n\t    for (var i = 0, n = arcs.length; i < n; ++i) arc(arcs[i], points);\n\t    return points;\n\t  }\n\n\t  function polygon(arcs) {\n\t    return arcs.map(line);\n\t  }\n\n\t  function geometry(o) {\n\t    o = Object.create(o);\n\t    o.properties = o.properties; // TODO possible duplicate\n\t    o.coordinates = geometryType[o.type](o.arcs);\n\t    //type is in o's prototype, which will be lost by worker.postMessage\n\t    o.type = o.type\n\t    return o;\n\t  }\n\t  var geometryType = {\n\t    LineString: line,\n\t    MultiLineString: polygon,\n\t    Polygon: polygon,\n\t    MultiPolygon: function(arcs) { return arcs.map(polygon); }\n\t  };\n\n\t  return o.type === \"GeometryCollection\"\n\t  ? (o = Object.create(o), o.geometries = o.geometries.map(geometry), o)\n\t  : geometry(o);\n\t}\n\n\n\tfunction reverse(array, n) {\n\t  var t, j = array.length, i = j - n; while (i < --j) t = array[i], array[i++] = array[j], array[j] = t;\n\t}\n\n\tHelpers.properties = function(obj) {\n\t  return obj.properties || {};\n\t}\n\n\tHelpers.transformer = function(tf) {\n\t  var kx = tf.scale[0],\n\t  ky = tf.scale[1],\n\t  dx = tf.translate[0],\n\t  dy = tf.translate[1];\n\n\t  function transform(c) {\n\t    return [c[0] * kx + dx, c[1] * ky + dy];\n\t  }\n\n\t  transform.invert = function(c) {\n\t    return [(c[0] - dx) / kx, (c[1]- dy) / ky];\n\t  };\n\n\t  return transform;\n\t};\n\n\tmodule.exports = Helpers\n\n\n/***/ }\n/******/ ]);", __webpack_require__.p + "2623f394d64d470ff1d8.worker.js");
	};

/***/ },
/* 2 */
/***/ function(module, exports) {

	// http://stackoverflow.com/questions/10343913/how-to-create-a-web-worker-from-a-string

	var URL = window.URL || window.webkitURL;
	module.exports = function(content, url) {
		try {
			try {
				var blob;
				try { // BlobBuilder = Deprecated, but widely implemented
					var BlobBuilder = window.BlobBuilder || window.WebKitBlobBuilder || window.MozBlobBuilder || window.MSBlobBuilder;
					blob = new BlobBuilder();
					blob.append(content);
					blob = blob.getBlob();
				} catch(e) { // The proposed API
					blob = new Blob([content]);
				}
				return new Worker(URL.createObjectURL(blob));
			} catch(e) {
				return new Worker('data:application/javascript,' + encodeURIComponent(content));
			}
		} catch(e) {
			return new Worker(url);
		}
	}

/***/ },
/* 3 */
/***/ function(module, exports, __webpack_require__) {

	//Use a partial d3 build, that doesn't need the DOM.
	var d3f = __webpack_require__(4)
	var Helpers = __webpack_require__(5)


	if (typeof onmessage !== 'undefined'){
	  onmessage = messaged
	} else {
	  function Wo(){}
	  Wo.prototype.postMessage = messaged
	  Wo.prototype.terminate = function(){}
	  module.exports = Wo
	}

	function messaged(event) {
	  var data = event
	  if (typeof event.data!== 'undefined') data = event.data

	  if (data.do === 'carto'){

	    var geo = data.geo,
	          topology = geo.topology,
	          geometries = geo.geometries,
	          path = geo.path,
	          translation = geo.translation,
	          projectionName = geo.projection.name,
	          scaling = geo.projection.scaling,
	          center = geo.projection.center,
	          translation = geo.projection.translation;

	    var values = data.values,
	        featureProperty = data.featureProperty,
	        task = data.task;



	    // copy it first
	    topology = Helpers.copy(topology);



	    // objects are projected into screen coordinates
	    // project the arcs into screen space


	    var projection =
	          d3f.geo[projectionName]()
	            .scale(scaling);

	    if (center != null) projection.center(center)
	    if (translation != null) projection.translate(translation)

	    var tf = Helpers.transformer(topology.transform),x,y,nArcVertices,vI,out1,nArcs=topology.arcs.length,aI=0,
	    projectedArcs = new Array(nArcs);
	    while(aI < nArcs){
	      x = 0;
	      y = 0;
	      nArcVertices = topology.arcs[aI].length;
	      vI = 0;
	      out1 = new Array(nArcVertices);
	      while( vI < nArcVertices){
	        topology.arcs[aI][vI][0] = (x += topology.arcs[aI][vI][0]);
	        topology.arcs[aI][vI][1] = (y += topology.arcs[aI][vI][1]);
	        out1[vI] = projection(tf(topology.arcs[aI][vI]));
	        vI++;
	      }
	      projectedArcs[aI++]=out1;

	    }

	    // path with identity projection
	    var path = d3f.geo.path()
	    .projection(null);


	    var objects = Helpers.object(projectedArcs, {type: "GeometryCollection", geometries: geometries})
	    .geometries.map(function(geom) {
	      return {
	        type: "Feature",
	        id: geom.id,
	        properties: Helpers.properties.call(null, geom, topology),
	        geometry: geom
	      };
	    });

	    function value(d){
	      return values[d.properties[featureProperty]]
	    }

	    var objectValues = objects.map(value),
	      totalValue = objectValues.reduce(function(a,b){return a + b;});

	    var iterations = 8;


	    //console.time("processing:" + task)
	    var i = 0;
	    while (i++ < iterations) {

	      //var areas = objects.map(path.area)
	      //var totalArea = areas.reduce(function(a,b){return a + b}),
	      var areas = [], totalArea = 0;
	      for (var k = 0; k < objects.length; k++){
	        var area = path.area(objects[k])
	        areas.push(area)
	        totalArea += area
	      }

	      var sizeErrorsTot = 0,
	      sizeErrorsNum = 0;

	      ///for i = 1 to n do
	      var meta = []
	      for (var j = 0; j < objects.length; j++){
	        var o = objects[j],
	            area = Math.abs(areas[j]), // XXX: why do we have negative areas?
	            v = +objectValues[j],
	            ///Compute AD i , the desired area of the ith cell
	            desired = totalArea * v / totalValue,
	            radius = Math.sqrt(area / Math.PI),
	            mass = Math.sqrt(desired / Math.PI) - radius,
	            sizeError = Math.max(area, desired) / Math.min(area, desired);

	        sizeErrorsTot+=sizeError;
	        sizeErrorsNum++;
	        // console.log(o.id, "@", j, "area:", area, "value:", v, "->", desired, radius, mass, sizeError);
	        meta.push({
	          id:         o.id,
	          area:       area,
	          centroid:   path.centroid(o),
	          value:      v,
	          desired:    desired,
	          range: 100 * (Math.abs(desired - area)) / (Math.sqrt(Math.PI * area)),
	          radius:     radius,
	          mass:       mass,
	          sizeError:  sizeError
	        })
	      }

	      var sizeError = sizeErrorsTot / sizeErrorsNum,
	          forceReductionFactor = 1 / (1 + sizeError);

	      // console.log("meta:", meta);
	      // console.log("  total area:", totalArea);
	      // console.log("  force reduction factor:", forceReductionFactor, "mean error:", sizeError);

	      var nArcVertices,vI,delta,nArcs=projectedArcs.length,aI=0,delta,nPolygon,pI,centroid,mass,radius,rSquared,dx,dy,distSquared,dist,Fij;
	      ///For each boundary line
	      while(aI < nArcs){
	        nArcVertices=projectedArcs[aI].length;
	        vI=0;
	        ///For each coordinate pair
	        while(vI < nArcVertices){
	          // create an array of vectors: [x, y]
	          delta = [0,0];
	          nPolygon = meta.length;
	          pI=0;
	          ///For each polygon centroid
	          while(pI < nPolygon) {
	            centroid =  meta[pI].centroid;
	            mass =      meta[pI].mass;
	            radius =    meta[pI].radius;
	            rSquared = (radius*radius);
	            dx = projectedArcs[aI][vI][0] - centroid[0];
	            dy = projectedArcs[aI][vI][1] - centroid[1];
	            distSquared = dx * dx + dy * dy;
	            dist=Math.sqrt(distSquared);
	            if (dist < meta[pI].range){
	              Fij = (dist > radius)
	              ? mass * radius / dist
	              : mass *
	              (distSquared / rSquared) *
	              (4 - 3 * dist / radius);
	              var tans = Helpers.arctans(dy, dx)
	              delta[0]+=(Fij * tans.cos);
	              delta[1]+=(Fij * tans.sin);
	            }
	            pI++;
	          }
	          projectedArcs[aI][vI][0] += (delta[0]*forceReductionFactor);
	          projectedArcs[aI][vI][1] += (delta[1]*forceReductionFactor);
	          vI++;
	        }
	        aI++;
	      }

	      // break if we hit the target size error
	      if (sizeError <= 1) break;
	    }

	    //console.timeEnd("processing:" + task)

	    var response = {
	      done: 'processing',
	      //geohson featureCollection
	      features: objects,
	      //arcs can be useful to reconstruct topojson :
	      //arcs: projectedArcs,
	      task: task
	    }

	    if (typeof self !== 'undefined'){
	      self.postMessage(response)
	    } else {
	      this.onmessage(response)
	    }
	  }
	}


/***/ },
/* 4 */
/***/ function(module, exports, __webpack_require__) {

	var __WEBPACK_AMD_DEFINE_FACTORY__, __WEBPACK_AMD_DEFINE_RESULT__;!function(){
	  var d3f = {version: "3.5.5"}; // semver
	function d3f_identity(d) {
	  return d;
	}
	var ε = 1e-6,
	    ε2 = ε * ε,
	    π = Math.PI,
	    τ = 2 * π,
	    τε = τ - ε,
	    halfπ = π / 2,
	    d3f_radians = π / 180,
	    d3f_degrees = 180 / π;

	function d3f_sgn(x) {
	  return x > 0 ? 1 : x < 0 ? -1 : 0;
	}

	// Returns the 2D cross product of AB and AC vectors, i.e., the z-component of
	// the 3D cross product in a quadrant I Cartesian coordinate system (+x is
	// right, +y is up). Returns a positive value if ABC is counter-clockwise,
	// negative if clockwise, and zero if the points are collinear.
	function d3f_cross2d(a, b, c) {
	  return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
	}

	function d3f_acos(x) {
	  return x > 1 ? 0 : x < -1 ? π : Math.acos(x);
	}

	function d3f_asin(x) {
	  return x > 1 ? halfπ : x < -1 ? -halfπ : Math.asin(x);
	}

	function d3f_sinh(x) {
	  return ((x = Math.exp(x)) - 1 / x) / 2;
	}

	function d3f_cosh(x) {
	  return ((x = Math.exp(x)) + 1 / x) / 2;
	}

	function d3f_tanh(x) {
	  return ((x = Math.exp(2 * x)) - 1) / (x + 1);
	}

	function d3f_haversin(x) {
	  return (x = Math.sin(x / 2)) * x;
	}
	d3f.geo = {};
	// Copies a variable number of methods from source to target.
	d3f.rebind = function(target, source) {
	  var i = 1, n = arguments.length, method;
	  while (++i < n) target[method = arguments[i]] = d3f_rebind(target, source, source[method]);
	  return target;
	};

	// Method is assumed to be a standard D3 getter-setter:
	// If passed with no arguments, gets the value.
	// If passed with arguments, sets the value and returns the target.
	function d3f_rebind(target, source, method) {
	  return function() {
	    var value = method.apply(source, arguments);
	    return value === source ? target : value;
	  };
	}
	function d3f_true() {
	  return true;
	}
	var abs = Math.abs;
	d3f.merge = function(arrays) {
	  var n = arrays.length,
	      m,
	      i = -1,
	      j = 0,
	      merged,
	      array;

	  while (++i < n) j += arrays[i].length;
	  merged = new Array(j);

	  while (--n >= 0) {
	    array = arrays[n];
	    m = array.length;
	    while (--m >= 0) {
	      merged[--j] = array[m];
	    }
	  }

	  return merged;
	};
	function d3f_noop() {}

	function d3f_geo_spherical(cartesian) {
	  return [
	    Math.atan2(cartesian[1], cartesian[0]),
	    d3f_asin(cartesian[2])
	  ];
	}

	function d3f_geo_sphericalEqual(a, b) {
	  return abs(a[0] - b[0]) < ε && abs(a[1] - b[1]) < ε;
	}

	// General spherical polygon clipping algorithm: takes a polygon, cuts it into
	// visible line segments and rejoins the segments by interpolating along the
	// clip edge.
	function d3f_geo_clipPolygon(segments, compare, clipStartInside, interpolate, listener) {
	  var subject = [],
	      clip = [];

	  segments.forEach(function(segment) {
	    if ((n = segment.length - 1) <= 0) return;
	    var n, p0 = segment[0], p1 = segment[n];

	    // If the first and last points of a segment are coincident, then treat as
	    // a closed ring.
	    // TODO if all rings are closed, then the winding order of the exterior
	    // ring should be checked.
	    if (d3f_geo_sphericalEqual(p0, p1)) {
	      listener.lineStart();
	      for (var i = 0; i < n; ++i) listener.point((p0 = segment[i])[0], p0[1]);
	      listener.lineEnd();
	      return;
	    }

	    var a = new d3f_geo_clipPolygonIntersection(p0, segment, null, true),
	        b = new d3f_geo_clipPolygonIntersection(p0, null, a, false);
	    a.o = b;
	    subject.push(a);
	    clip.push(b);
	    a = new d3f_geo_clipPolygonIntersection(p1, segment, null, false);
	    b = new d3f_geo_clipPolygonIntersection(p1, null, a, true);
	    a.o = b;
	    subject.push(a);
	    clip.push(b);
	  });
	  clip.sort(compare);
	  d3f_geo_clipPolygonLinkCircular(subject);
	  d3f_geo_clipPolygonLinkCircular(clip);
	  if (!subject.length) return;

	  for (var i = 0, entry = clipStartInside, n = clip.length; i < n; ++i) {
	    clip[i].e = entry = !entry;
	  }

	  var start = subject[0],
	      points,
	      point;
	  while (1) {
	    // Find first unvisited intersection.
	    var current = start,
	        isSubject = true;
	    while (current.v) if ((current = current.n) === start) return;
	    points = current.z;
	    listener.lineStart();
	    do {
	      current.v = current.o.v = true;
	      if (current.e) {
	        if (isSubject) {
	          for (var i = 0, n = points.length; i < n; ++i) listener.point((point = points[i])[0], point[1]);
	        } else {
	          interpolate(current.x, current.n.x, 1, listener);
	        }
	        current = current.n;
	      } else {
	        if (isSubject) {
	          points = current.p.z;
	          for (var i = points.length - 1; i >= 0; --i) listener.point((point = points[i])[0], point[1]);
	        } else {
	          interpolate(current.x, current.p.x, -1, listener);
	        }
	        current = current.p;
	      }
	      current = current.o;
	      points = current.z;
	      isSubject = !isSubject;
	    } while (!current.v);
	    listener.lineEnd();
	  }
	}

	function d3f_geo_clipPolygonLinkCircular(array) {
	  if (!(n = array.length)) return;
	  var n,
	      i = 0,
	      a = array[0],
	      b;
	  while (++i < n) {
	    a.n = b = array[i];
	    b.p = a;
	    a = b;
	  }
	  a.n = b = array[0];
	  b.p = a;
	}

	function d3f_geo_clipPolygonIntersection(point, points, other, entry) {
	  this.x = point;
	  this.z = points;
	  this.o = other; // another intersection
	  this.e = entry; // is an entry?
	  this.v = false; // visited
	  this.n = this.p = null; // next & previous
	}

	function d3f_geo_clip(pointVisible, clipLine, interpolate, clipStart) {
	  return function(rotate, listener) {
	    var line = clipLine(listener),
	        rotatedClipStart = rotate.invert(clipStart[0], clipStart[1]);

	    var clip = {
	      point: point,
	      lineStart: lineStart,
	      lineEnd: lineEnd,
	      polygonStart: function() {
	        clip.point = pointRing;
	        clip.lineStart = ringStart;
	        clip.lineEnd = ringEnd;
	        segments = [];
	        polygon = [];
	      },
	      polygonEnd: function() {
	        clip.point = point;
	        clip.lineStart = lineStart;
	        clip.lineEnd = lineEnd;

	        segments = d3f.merge(segments);
	        var clipStartInside = d3f_geo_pointInPolygon(rotatedClipStart, polygon);
	        if (segments.length) {
	          if (!polygonStarted) listener.polygonStart(), polygonStarted = true;
	          d3f_geo_clipPolygon(segments, d3f_geo_clipSort, clipStartInside, interpolate, listener);
	        } else if (clipStartInside) {
	          if (!polygonStarted) listener.polygonStart(), polygonStarted = true;
	          listener.lineStart();
	          interpolate(null, null, 1, listener);
	          listener.lineEnd();
	        }
	        if (polygonStarted) listener.polygonEnd(), polygonStarted = false;
	        segments = polygon = null;
	      },
	      sphere: function() {
	        listener.polygonStart();
	        listener.lineStart();
	        interpolate(null, null, 1, listener);
	        listener.lineEnd();
	        listener.polygonEnd();
	      }
	    };

	    function point(λ, φ) {
	      var point = rotate(λ, φ);
	      if (pointVisible(λ = point[0], φ = point[1])) listener.point(λ, φ);
	    }
	    function pointLine(λ, φ) {
	      var point = rotate(λ, φ);
	      line.point(point[0], point[1]);
	    }
	    function lineStart() { clip.point = pointLine; line.lineStart(); }
	    function lineEnd() { clip.point = point; line.lineEnd(); }

	    var segments;

	    var buffer = d3f_geo_clipBufferListener(),
	        ringListener = clipLine(buffer),
	        polygonStarted = false,
	        polygon,
	        ring;

	    function pointRing(λ, φ) {
	      ring.push([λ, φ]);
	      var point = rotate(λ, φ);
	      ringListener.point(point[0], point[1]);
	    }

	    function ringStart() {
	      ringListener.lineStart();
	      ring = [];
	    }

	    function ringEnd() {
	      pointRing(ring[0][0], ring[0][1]);
	      ringListener.lineEnd();

	      var clean = ringListener.clean(),
	          ringSegments = buffer.buffer(),
	          segment,
	          n = ringSegments.length;

	      ring.pop();
	      polygon.push(ring);
	      ring = null;

	      if (!n) return;

	      // No intersections.
	      if (clean & 1) {
	        segment = ringSegments[0];
	        var n = segment.length - 1,
	            i = -1,
	            point;
	        if (n > 0) {
	          if (!polygonStarted) listener.polygonStart(), polygonStarted = true;
	          listener.lineStart();
	          while (++i < n) listener.point((point = segment[i])[0], point[1]);
	          listener.lineEnd();
	        }
	        return;
	      }

	      // Rejoin connected segments.
	      // TODO reuse bufferListener.rejoin()?
	      if (n > 1 && clean & 2) ringSegments.push(ringSegments.pop().concat(ringSegments.shift()));

	      segments.push(ringSegments.filter(d3f_geo_clipSegmentLength1));
	    }

	    return clip;
	  };
	}

	function d3f_geo_clipSegmentLength1(segment) {
	  return segment.length > 1;
	}

	function d3f_geo_clipBufferListener() {
	  var lines = [],
	      line;
	  return {
	    lineStart: function() { lines.push(line = []); },
	    point: function(λ, φ) { line.push([λ, φ]); },
	    lineEnd: d3f_noop,
	    buffer: function() {
	      var buffer = lines;
	      lines = [];
	      line = null;
	      return buffer;
	    },
	    rejoin: function() {
	      if (lines.length > 1) lines.push(lines.pop().concat(lines.shift()));
	    }
	  };
	}

	// Intersection points are sorted along the clip edge. For both antimeridian
	// cutting and circle clipping, the same comparison is used.
	function d3f_geo_clipSort(a, b) {
	  return ((a = a.x)[0] < 0 ? a[1] - halfπ - ε : halfπ - a[1])
	       - ((b = b.x)[0] < 0 ? b[1] - halfπ - ε : halfπ - b[1]);
	}

	var d3f_geo_clipAntimeridian = d3f_geo_clip(
	    d3f_true,
	    d3f_geo_clipAntimeridianLine,
	    d3f_geo_clipAntimeridianInterpolate,
	    [-π, -π / 2]);

	// Takes a line and cuts into visible segments. Return values:
	//   0: there were intersections or the line was empty.
	//   1: no intersections.
	//   2: there were intersections, and the first and last segments should be
	//      rejoined.
	function d3f_geo_clipAntimeridianLine(listener) {
	  var λ0 = NaN,
	      φ0 = NaN,
	      sλ0 = NaN,
	      clean; // no intersections

	  return {
	    lineStart: function() {
	      listener.lineStart();
	      clean = 1;
	    },
	    point: function(λ1, φ1) {
	      var sλ1 = λ1 > 0 ? π : -π,
	          dλ = abs(λ1 - λ0);
	      if (abs(dλ - π) < ε) { // line crosses a pole
	        listener.point(λ0, φ0 = (φ0 + φ1) / 2 > 0 ? halfπ : -halfπ);
	        listener.point(sλ0, φ0);
	        listener.lineEnd();
	        listener.lineStart();
	        listener.point(sλ1, φ0);
	        listener.point(λ1, φ0);
	        clean = 0;
	      } else if (sλ0 !== sλ1 && dλ >= π) { // line crosses antimeridian
	        // handle degeneracies
	        if (abs(λ0 - sλ0) < ε) λ0 -= sλ0 * ε;
	        if (abs(λ1 - sλ1) < ε) λ1 -= sλ1 * ε;
	        φ0 = d3f_geo_clipAntimeridianIntersect(λ0, φ0, λ1, φ1);
	        listener.point(sλ0, φ0);
	        listener.lineEnd();
	        listener.lineStart();
	        listener.point(sλ1, φ0);
	        clean = 0;
	      }
	      listener.point(λ0 = λ1, φ0 = φ1);
	      sλ0 = sλ1;
	    },
	    lineEnd: function() {
	      listener.lineEnd();
	      λ0 = φ0 = NaN;
	    },
	    // if there are intersections, we always rejoin the first and last segments.
	    clean: function() { return 2 - clean; }
	  };
	}

	function d3f_geo_clipAntimeridianIntersect(λ0, φ0, λ1, φ1) {
	  var cosφ0,
	      cosφ1,
	      sinλ0_λ1 = Math.sin(λ0 - λ1);
	  return abs(sinλ0_λ1) > ε
	      ? Math.atan((Math.sin(φ0) * (cosφ1 = Math.cos(φ1)) * Math.sin(λ1)
	                 - Math.sin(φ1) * (cosφ0 = Math.cos(φ0)) * Math.sin(λ0))
	                 / (cosφ0 * cosφ1 * sinλ0_λ1))
	      : (φ0 + φ1) / 2;
	}

	function d3f_geo_clipAntimeridianInterpolate(from, to, direction, listener) {
	  var φ;
	  if (from == null) {
	    φ = direction * halfπ;
	    listener.point(-π,  φ);
	    listener.point( 0,  φ);
	    listener.point( π,  φ);
	    listener.point( π,  0);
	    listener.point( π, -φ);
	    listener.point( 0, -φ);
	    listener.point(-π, -φ);
	    listener.point(-π,  0);
	    listener.point(-π,  φ);
	  } else if (abs(from[0] - to[0]) > ε) {
	    var s = from[0] < to[0] ? π : -π;
	    φ = direction * s / 2;
	    listener.point(-s, φ);
	    listener.point( 0, φ);
	    listener.point( s, φ);
	  } else {
	    listener.point(to[0], to[1]);
	  }
	}
	// TODO
	// cross and scale return new vectors,
	// whereas add and normalize operate in-place

	function d3f_geo_cartesian(spherical) {
	  var λ = spherical[0],
	      φ = spherical[1],
	      cosφ = Math.cos(φ);
	  return [
	    cosφ * Math.cos(λ),
	    cosφ * Math.sin(λ),
	    Math.sin(φ)
	  ];
	}

	function d3f_geo_cartesianDot(a, b) {
	  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}

	function d3f_geo_cartesianCross(a, b) {
	  return [
	    a[1] * b[2] - a[2] * b[1],
	    a[2] * b[0] - a[0] * b[2],
	    a[0] * b[1] - a[1] * b[0]
	  ];
	}

	function d3f_geo_cartesianAdd(a, b) {
	  a[0] += b[0];
	  a[1] += b[1];
	  a[2] += b[2];
	}

	function d3f_geo_cartesianScale(vector, k) {
	  return [
	    vector[0] * k,
	    vector[1] * k,
	    vector[2] * k
	  ];
	}

	function d3f_geo_cartesianNormalize(d) {
	  var l = Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
	  d[0] /= l;
	  d[1] /= l;
	  d[2] /= l;
	}
	function d3f_geo_compose(a, b) {

	  function compose(x, y) {
	    return x = a(x, y), b(x[0], x[1]);
	  }

	  if (a.invert && b.invert) compose.invert = function(x, y) {
	    return x = b.invert(x, y), x && a.invert(x[0], x[1]);
	  };

	  return compose;
	}

	function d3f_geo_equirectangular(λ, φ) {
	  return [λ, φ];
	}

	(d3f.geo.equirectangular = function() {
	  return d3f_geo_projection(d3f_geo_equirectangular);
	}).raw = d3f_geo_equirectangular.invert = d3f_geo_equirectangular;

	d3f.geo.rotation = function(rotate) {
	  rotate = d3f_geo_rotation(rotate[0] % 360 * d3f_radians, rotate[1] * d3f_radians, rotate.length > 2 ? rotate[2] * d3f_radians : 0);

	  function forward(coordinates) {
	    coordinates = rotate(coordinates[0] * d3f_radians, coordinates[1] * d3f_radians);
	    return coordinates[0] *= d3f_degrees, coordinates[1] *= d3f_degrees, coordinates;
	  }

	  forward.invert = function(coordinates) {
	    coordinates = rotate.invert(coordinates[0] * d3f_radians, coordinates[1] * d3f_radians);
	    return coordinates[0] *= d3f_degrees, coordinates[1] *= d3f_degrees, coordinates;
	  };

	  return forward;
	};

	function d3f_geo_identityRotation(λ, φ) {
	  return [λ > π ? λ - τ : λ < -π ? λ + τ : λ, φ];
	}

	d3f_geo_identityRotation.invert = d3f_geo_equirectangular;

	// Note: |δλ| must be < 2π
	function d3f_geo_rotation(δλ, δφ, δγ) {
	  return δλ ? (δφ || δγ ? d3f_geo_compose(d3f_geo_rotationλ(δλ), d3f_geo_rotationφγ(δφ, δγ))
	    : d3f_geo_rotationλ(δλ))
	    : (δφ || δγ ? d3f_geo_rotationφγ(δφ, δγ)
	    : d3f_geo_identityRotation);
	}

	function d3f_geo_forwardRotationλ(δλ) {
	  return function(λ, φ) {
	    return λ += δλ, [λ > π ? λ - τ : λ < -π ? λ + τ : λ, φ];
	  };
	}

	function d3f_geo_rotationλ(δλ) {
	  var rotation = d3f_geo_forwardRotationλ(δλ);
	  rotation.invert = d3f_geo_forwardRotationλ(-δλ);
	  return rotation;
	}

	function d3f_geo_rotationφγ(δφ, δγ) {
	  var cosδφ = Math.cos(δφ),
	      sinδφ = Math.sin(δφ),
	      cosδγ = Math.cos(δγ),
	      sinδγ = Math.sin(δγ);

	  function rotation(λ, φ) {
	    var cosφ = Math.cos(φ),
	        x = Math.cos(λ) * cosφ,
	        y = Math.sin(λ) * cosφ,
	        z = Math.sin(φ),
	        k = z * cosδφ + x * sinδφ;
	    return [
	      Math.atan2(y * cosδγ - k * sinδγ, x * cosδφ - z * sinδφ),
	      d3f_asin(k * cosδγ + y * sinδγ)
	    ];
	  }

	  rotation.invert = function(λ, φ) {
	    var cosφ = Math.cos(φ),
	        x = Math.cos(λ) * cosφ,
	        y = Math.sin(λ) * cosφ,
	        z = Math.sin(φ),
	        k = z * cosδγ - y * sinδγ;
	    return [
	      Math.atan2(y * cosδγ + z * sinδγ, x * cosδφ + k * sinδφ),
	      d3f_asin(k * cosδφ - x * sinδφ)
	    ];
	  };

	  return rotation;
	}

	d3f.geo.circle = function() {
	  var origin = [0, 0],
	      angle,
	      precision = 6,
	      interpolate;

	  function circle() {
	    var center = typeof origin === "function" ? origin.apply(this, arguments) : origin,
	        rotate = d3f_geo_rotation(-center[0] * d3f_radians, -center[1] * d3f_radians, 0).invert,
	        ring = [];

	    interpolate(null, null, 1, {
	      point: function(x, y) {
	        ring.push(x = rotate(x, y));
	        x[0] *= d3f_degrees, x[1] *= d3f_degrees;
	      }
	    });

	    return {type: "Polygon", coordinates: [ring]};
	  }

	  circle.origin = function(x) {
	    if (!arguments.length) return origin;
	    origin = x;
	    return circle;
	  };

	  circle.angle = function(x) {
	    if (!arguments.length) return angle;
	    interpolate = d3f_geo_circleInterpolate((angle = +x) * d3f_radians, precision * d3f_radians);
	    return circle;
	  };

	  circle.precision = function(_) {
	    if (!arguments.length) return precision;
	    interpolate = d3f_geo_circleInterpolate(angle * d3f_radians, (precision = +_) * d3f_radians);
	    return circle;
	  };

	  return circle.angle(90);
	};

	// Interpolates along a circle centered at [0°, 0°], with a given radius and
	// precision.
	function d3f_geo_circleInterpolate(radius, precision) {
	  var cr = Math.cos(radius),
	      sr = Math.sin(radius);
	  return function(from, to, direction, listener) {
	    var step = direction * precision;
	    if (from != null) {
	      from = d3f_geo_circleAngle(cr, from);
	      to = d3f_geo_circleAngle(cr, to);
	      if (direction > 0 ? from < to: from > to) from += direction * τ;
	    } else {
	      from = radius + direction * τ;
	      to = radius - .5 * step;
	    }
	    for (var point, t = from; direction > 0 ? t > to : t < to; t -= step) {
	      listener.point((point = d3f_geo_spherical([
	        cr,
	        -sr * Math.cos(t),
	        -sr * Math.sin(t)
	      ]))[0], point[1]);
	    }
	  };
	}

	// Signed angle of a cartesian point relative to [cr, 0, 0].
	function d3f_geo_circleAngle(cr, point) {
	  var a = d3f_geo_cartesian(point);
	  a[0] -= cr;
	  d3f_geo_cartesianNormalize(a);
	  var angle = d3f_acos(-a[1]);
	  return ((-a[2] < 0 ? -angle : angle) + 2 * Math.PI - ε) % (2 * Math.PI);
	}
	// Adds floating point numbers with twice the normal precision.
	// Reference: J. R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and
	// Fast Robust Geometric Predicates, Discrete & Computational Geometry 18(3)
	// 305–363 (1997).
	// Code adapted from GeographicLib by Charles F. F. Karney,
	// http://geographiclib.sourceforge.net/
	// See lib/geographiclib/LICENSE for details.

	function d3f_adder() {}

	d3f_adder.prototype = {
	  s: 0, // rounded value
	  t: 0, // exact error
	  add: function(y) {
	    d3f_adderSum(y, this.t, d3f_adderTemp);
	    d3f_adderSum(d3f_adderTemp.s, this.s, this);
	    if (this.s) this.t += d3f_adderTemp.t;
	    else this.s = d3f_adderTemp.t;
	  },
	  reset: function() {
	    this.s = this.t = 0;
	  },
	  valueOf: function() {
	    return this.s;
	  }
	};

	var d3f_adderTemp = new d3f_adder;

	function d3f_adderSum(a, b, o) {
	  var x = o.s = a + b, // a + b
	      bv = x - a, av = x - bv; // b_virtual & a_virtual
	  o.t = (a - av) + (b - bv); // a_roundoff + b_roundoff
	}

	d3f.geo.stream = function(object, listener) {
	  if (object && d3f_geo_streamObjectType.hasOwnProperty(object.type)) {
	    d3f_geo_streamObjectType[object.type](object, listener);
	  } else {
	    d3f_geo_streamGeometry(object, listener);
	  }
	};

	function d3f_geo_streamGeometry(geometry, listener) {
	  if (geometry && d3f_geo_streamGeometryType.hasOwnProperty(geometry.type)) {
	    d3f_geo_streamGeometryType[geometry.type](geometry, listener);
	  }
	}

	var d3f_geo_streamObjectType = {
	  Feature: function(feature, listener) {
	    d3f_geo_streamGeometry(feature.geometry, listener);
	  },
	  FeatureCollection: function(object, listener) {
	    var features = object.features, i = -1, n = features.length;
	    while (++i < n) d3f_geo_streamGeometry(features[i].geometry, listener);
	  }
	};

	var d3f_geo_streamGeometryType = {
	  Sphere: function(object, listener) {
	    listener.sphere();
	  },
	  Point: function(object, listener) {
	    object = object.coordinates;
	    listener.point(object[0], object[1], object[2]);
	  },
	  MultiPoint: function(object, listener) {
	    var coordinates = object.coordinates, i = -1, n = coordinates.length;
	    while (++i < n) object = coordinates[i], listener.point(object[0], object[1], object[2]);
	  },
	  LineString: function(object, listener) {
	    d3f_geo_streamLine(object.coordinates, listener, 0);
	  },
	  MultiLineString: function(object, listener) {
	    var coordinates = object.coordinates, i = -1, n = coordinates.length;
	    while (++i < n) d3f_geo_streamLine(coordinates[i], listener, 0);
	  },
	  Polygon: function(object, listener) {
	    d3f_geo_streamPolygon(object.coordinates, listener);
	  },
	  MultiPolygon: function(object, listener) {
	    var coordinates = object.coordinates, i = -1, n = coordinates.length;
	    while (++i < n) d3f_geo_streamPolygon(coordinates[i], listener);
	  },
	  GeometryCollection: function(object, listener) {
	    var geometries = object.geometries, i = -1, n = geometries.length;
	    while (++i < n) d3f_geo_streamGeometry(geometries[i], listener);
	  }
	};

	function d3f_geo_streamLine(coordinates, listener, closed) {
	  var i = -1, n = coordinates.length - closed, coordinate;
	  listener.lineStart();
	  while (++i < n) coordinate = coordinates[i], listener.point(coordinate[0], coordinate[1], coordinate[2]);
	  listener.lineEnd();
	}

	function d3f_geo_streamPolygon(coordinates, listener) {
	  var i = -1, n = coordinates.length;
	  listener.polygonStart();
	  while (++i < n) d3f_geo_streamLine(coordinates[i], listener, 1);
	  listener.polygonEnd();
	}

	d3f.geo.area = function(object) {
	  d3f_geo_areaSum = 0;
	  d3f.geo.stream(object, d3f_geo_area);
	  return d3f_geo_areaSum;
	};

	var d3f_geo_areaSum,
	    d3f_geo_areaRingSum = new d3f_adder;

	var d3f_geo_area = {
	  sphere: function() { d3f_geo_areaSum += 4 * π; },
	  point: d3f_noop,
	  lineStart: d3f_noop,
	  lineEnd: d3f_noop,

	  // Only count area for polygon rings.
	  polygonStart: function() {
	    d3f_geo_areaRingSum.reset();
	    d3f_geo_area.lineStart = d3f_geo_areaRingStart;
	  },
	  polygonEnd: function() {
	    var area = 2 * d3f_geo_areaRingSum;
	    d3f_geo_areaSum += area < 0 ? 4 * π + area : area;
	    d3f_geo_area.lineStart = d3f_geo_area.lineEnd = d3f_geo_area.point = d3f_noop;
	  }
	};

	function d3f_geo_areaRingStart() {
	  var λ00, φ00, λ0, cosφ0, sinφ0; // start point and previous point

	  // For the first point, …
	  d3f_geo_area.point = function(λ, φ) {
	    d3f_geo_area.point = nextPoint;
	    λ0 = (λ00 = λ) * d3f_radians, cosφ0 = Math.cos(φ = (φ00 = φ) * d3f_radians / 2 + π / 4), sinφ0 = Math.sin(φ);
	  };

	  // For subsequent points, …
	  function nextPoint(λ, φ) {
	    λ *= d3f_radians;
	    φ = φ * d3f_radians / 2 + π / 4; // half the angular distance from south pole

	    // Spherical excess E for a spherical triangle with vertices: south pole,
	    // previous point, current point.  Uses a formula derived from Cagnoli’s
	    // theorem.  See Todhunter, Spherical Trig. (1871), Sec. 103, Eq. (2).
	    var dλ = λ - λ0,
	        sdλ = dλ >= 0 ? 1 : -1,
	        adλ = sdλ * dλ,
	        cosφ = Math.cos(φ),
	        sinφ = Math.sin(φ),
	        k = sinφ0 * sinφ,
	        u = cosφ0 * cosφ + k * Math.cos(adλ),
	        v = k * sdλ * Math.sin(adλ);
	    d3f_geo_areaRingSum.add(Math.atan2(v, u));

	    // Advance the previous points.
	    λ0 = λ, cosφ0 = cosφ, sinφ0 = sinφ;
	  }

	  // For the last point, return to the start.
	  d3f_geo_area.lineEnd = function() {
	    nextPoint(λ00, φ00);
	  };
	}

	function d3f_geo_pointInPolygon(point, polygon) {
	  var meridian = point[0],
	      parallel = point[1],
	      meridianNormal = [Math.sin(meridian), -Math.cos(meridian), 0],
	      polarAngle = 0,
	      winding = 0;
	  d3f_geo_areaRingSum.reset();

	  for (var i = 0, n = polygon.length; i < n; ++i) {
	    var ring = polygon[i],
	        m = ring.length;
	    if (!m) continue;
	    var point0 = ring[0],
	        λ0 = point0[0],
	        φ0 = point0[1] / 2 + π / 4,
	        sinφ0 = Math.sin(φ0),
	        cosφ0 = Math.cos(φ0),
	        j = 1;

	    while (true) {
	      if (j === m) j = 0;
	      point = ring[j];
	      var λ = point[0],
	          φ = point[1] / 2 + π / 4,
	          sinφ = Math.sin(φ),
	          cosφ = Math.cos(φ),
	          dλ = λ - λ0,
	          sdλ = dλ >= 0 ? 1 : -1,
	          adλ = sdλ * dλ,
	          antimeridian = adλ > π,
	          k = sinφ0 * sinφ;
	      d3f_geo_areaRingSum.add(Math.atan2(k * sdλ * Math.sin(adλ), cosφ0 * cosφ + k * Math.cos(adλ)));

	      polarAngle += antimeridian ? dλ + sdλ * τ : dλ;

	      // Are the longitudes either side of the point's meridian, and are the
	      // latitudes smaller than the parallel?
	      if (antimeridian ^ λ0 >= meridian ^ λ >= meridian) {
	        var arc = d3f_geo_cartesianCross(d3f_geo_cartesian(point0), d3f_geo_cartesian(point));
	        d3f_geo_cartesianNormalize(arc);
	        var intersection = d3f_geo_cartesianCross(meridianNormal, arc);
	        d3f_geo_cartesianNormalize(intersection);
	        var φarc = (antimeridian ^ dλ >= 0 ? -1 : 1) * d3f_asin(intersection[2]);
	        if (parallel > φarc || parallel === φarc && (arc[0] || arc[1])) {
	          winding += antimeridian ^ dλ >= 0 ? 1 : -1;
	        }
	      }
	      if (!j++) break;
	      λ0 = λ, sinφ0 = sinφ, cosφ0 = cosφ, point0 = point;
	    }
	  }

	  // First, determine whether the South pole is inside or outside:
	  //
	  // It is inside if:
	  // * the polygon winds around it in a clockwise direction.
	  // * the polygon does not (cumulatively) wind around it, but has a negative
	  //   (counter-clockwise) area.
	  //
	  // Second, count the (signed) number of times a segment crosses a meridian
	  // from the point to the South pole.  If it is zero, then the point is the
	  // same side as the South pole.

	  return (polarAngle < -ε || polarAngle < ε && d3f_geo_areaRingSum < 0) ^ (winding & 1);
	}

	// Clip features against a small circle centered at [0°, 0°].
	function d3f_geo_clipCircle(radius) {
	  var cr = Math.cos(radius),
	      smallRadius = cr > 0,
	      notHemisphere = abs(cr) > ε, // TODO optimise for this common case
	      interpolate = d3f_geo_circleInterpolate(radius, 6 * d3f_radians);

	  return d3f_geo_clip(visible, clipLine, interpolate, smallRadius ? [0, -radius] : [-π, radius - π]);

	  function visible(λ, φ) {
	    return Math.cos(λ) * Math.cos(φ) > cr;
	  }

	  // Takes a line and cuts into visible segments. Return values used for
	  // polygon clipping:
	  //   0: there were intersections or the line was empty.
	  //   1: no intersections.
	  //   2: there were intersections, and the first and last segments should be
	  //      rejoined.
	  function clipLine(listener) {
	    var point0, // previous point
	        c0, // code for previous point
	        v0, // visibility of previous point
	        v00, // visibility of first point
	        clean; // no intersections
	    return {
	      lineStart: function() {
	        v00 = v0 = false;
	        clean = 1;
	      },
	      point: function(λ, φ) {
	        var point1 = [λ, φ],
	            point2,
	            v = visible(λ, φ),
	            c = smallRadius
	              ? v ? 0 : code(λ, φ)
	              : v ? code(λ + (λ < 0 ? π : -π), φ) : 0;
	        if (!point0 && (v00 = v0 = v)) listener.lineStart();
	        // Handle degeneracies.
	        // TODO ignore if not clipping polygons.
	        if (v !== v0) {
	          point2 = intersect(point0, point1);
	          if (d3f_geo_sphericalEqual(point0, point2) || d3f_geo_sphericalEqual(point1, point2)) {
	            point1[0] += ε;
	            point1[1] += ε;
	            v = visible(point1[0], point1[1]);
	          }
	        }
	        if (v !== v0) {
	          clean = 0;
	          if (v) {
	            // outside going in
	            listener.lineStart();
	            point2 = intersect(point1, point0);
	            listener.point(point2[0], point2[1]);
	          } else {
	            // inside going out
	            point2 = intersect(point0, point1);
	            listener.point(point2[0], point2[1]);
	            listener.lineEnd();
	          }
	          point0 = point2;
	        } else if (notHemisphere && point0 && smallRadius ^ v) {
	          var t;
	          // If the codes for two points are different, or are both zero,
	          // and there this segment intersects with the small circle.
	          if (!(c & c0) && (t = intersect(point1, point0, true))) {
	            clean = 0;
	            if (smallRadius) {
	              listener.lineStart();
	              listener.point(t[0][0], t[0][1]);
	              listener.point(t[1][0], t[1][1]);
	              listener.lineEnd();
	            } else {
	              listener.point(t[1][0], t[1][1]);
	              listener.lineEnd();
	              listener.lineStart();
	              listener.point(t[0][0], t[0][1]);
	            }
	          }
	        }
	        if (v && (!point0 || !d3f_geo_sphericalEqual(point0, point1))) {
	          listener.point(point1[0], point1[1]);
	        }
	        point0 = point1, v0 = v, c0 = c;
	      },
	      lineEnd: function() {
	        if (v0) listener.lineEnd();
	        point0 = null;
	      },
	      // Rejoin first and last segments if there were intersections and the first
	      // and last points were visible.
	      clean: function() { return clean | ((v00 && v0) << 1); }
	    };
	  }

	  // Intersects the great circle between a and b with the clip circle.
	  function intersect(a, b, two) {
	    var pa = d3f_geo_cartesian(a),
	        pb = d3f_geo_cartesian(b);

	    // We have two planes, n1.p = d1 and n2.p = d2.
	    // Find intersection line p(t) = c1 n1 + c2 n2 + t (n1 ⨯ n2).
	    var n1 = [1, 0, 0], // normal
	        n2 = d3f_geo_cartesianCross(pa, pb),
	        n2n2 = d3f_geo_cartesianDot(n2, n2),
	        n1n2 = n2[0], // d3f_geo_cartesianDot(n1, n2),
	        determinant = n2n2 - n1n2 * n1n2;

	    // Two polar points.
	    if (!determinant) return !two && a;

	    var c1 =  cr * n2n2 / determinant,
	        c2 = -cr * n1n2 / determinant,
	        n1xn2 = d3f_geo_cartesianCross(n1, n2),
	        A = d3f_geo_cartesianScale(n1, c1),
	        B = d3f_geo_cartesianScale(n2, c2);
	    d3f_geo_cartesianAdd(A, B);

	    // Solve |p(t)|^2 = 1.
	    var u = n1xn2,
	        w = d3f_geo_cartesianDot(A, u),
	        uu = d3f_geo_cartesianDot(u, u),
	        t2 = w * w - uu * (d3f_geo_cartesianDot(A, A) - 1);

	    if (t2 < 0) return;

	    var t = Math.sqrt(t2),
	        q = d3f_geo_cartesianScale(u, (-w - t) / uu);
	    d3f_geo_cartesianAdd(q, A);
	    q = d3f_geo_spherical(q);
	    if (!two) return q;

	    // Two intersection points.
	    var λ0 = a[0],
	        λ1 = b[0],
	        φ0 = a[1],
	        φ1 = b[1],
	        z;
	    if (λ1 < λ0) z = λ0, λ0 = λ1, λ1 = z;
	    var δλ = λ1 - λ0,
	        polar = abs(δλ - π) < ε,
	        meridian = polar || δλ < ε;

	    if (!polar && φ1 < φ0) z = φ0, φ0 = φ1, φ1 = z;

	    // Check that the first point is between a and b.
	    if (meridian
	        ? polar
	          ? φ0 + φ1 > 0 ^ q[1] < (abs(q[0] - λ0) < ε ? φ0 : φ1)
	          : φ0 <= q[1] && q[1] <= φ1
	        : δλ > π ^ (λ0 <= q[0] && q[0] <= λ1)) {
	      var q1 = d3f_geo_cartesianScale(u, (-w + t) / uu);
	      d3f_geo_cartesianAdd(q1, A);
	      return [q, d3f_geo_spherical(q1)];
	    }
	  }

	  // Generates a 4-bit vector representing the location of a point relative to
	  // the small circle's bounding box.
	  function code(λ, φ) {
	    var r = smallRadius ? radius : π - radius,
	        code = 0;
	    if (λ < -r) code |= 1; // left
	    else if (λ > r) code |= 2; // right
	    if (φ < -r) code |= 4; // below
	    else if (φ > r) code |= 8; // above
	    return code;
	  }
	}

	// Liang–Barsky line clipping.
	function d3f_geom_clipLine(x0, y0, x1, y1) {
	  return function(line) {
	    var a = line.a,
	        b = line.b,
	        ax = a.x,
	        ay = a.y,
	        bx = b.x,
	        by = b.y,
	        t0 = 0,
	        t1 = 1,
	        dx = bx - ax,
	        dy = by - ay,
	        r;

	    r = x0 - ax;
	    if (!dx && r > 0) return;
	    r /= dx;
	    if (dx < 0) {
	      if (r < t0) return;
	      if (r < t1) t1 = r;
	    } else if (dx > 0) {
	      if (r > t1) return;
	      if (r > t0) t0 = r;
	    }

	    r = x1 - ax;
	    if (!dx && r < 0) return;
	    r /= dx;
	    if (dx < 0) {
	      if (r > t1) return;
	      if (r > t0) t0 = r;
	    } else if (dx > 0) {
	      if (r < t0) return;
	      if (r < t1) t1 = r;
	    }

	    r = y0 - ay;
	    if (!dy && r > 0) return;
	    r /= dy;
	    if (dy < 0) {
	      if (r < t0) return;
	      if (r < t1) t1 = r;
	    } else if (dy > 0) {
	      if (r > t1) return;
	      if (r > t0) t0 = r;
	    }

	    r = y1 - ay;
	    if (!dy && r < 0) return;
	    r /= dy;
	    if (dy < 0) {
	      if (r > t1) return;
	      if (r > t0) t0 = r;
	    } else if (dy > 0) {
	      if (r < t0) return;
	      if (r < t1) t1 = r;
	    }

	    if (t0 > 0) line.a = {x: ax + t0 * dx, y: ay + t0 * dy};
	    if (t1 < 1) line.b = {x: ax + t1 * dx, y: ay + t1 * dy};
	    return line;
	  };
	}

	var d3f_geo_clipExtentMAX = 1e9;

	d3f.geo.clipExtent = function() {
	  var x0, y0, x1, y1,
	      stream,
	      clip,
	      clipExtent = {
	        stream: function(output) {
	          if (stream) stream.valid = false;
	          stream = clip(output);
	          stream.valid = true; // allow caching by d3f.geo.path
	          return stream;
	        },
	        extent: function(_) {
	          if (!arguments.length) return [[x0, y0], [x1, y1]];
	          clip = d3f_geo_clipExtent(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]);
	          if (stream) stream.valid = false, stream = null;
	          return clipExtent;
	        }
	      };
	  return clipExtent.extent([[0, 0], [960, 500]]);
	};

	function d3f_geo_clipExtent(x0, y0, x1, y1) {
	  return function(listener) {
	    var listener_ = listener,
	        bufferListener = d3f_geo_clipBufferListener(),
	        clipLine = d3f_geom_clipLine(x0, y0, x1, y1),
	        segments,
	        polygon,
	        ring;

	    var clip = {
	      point: point,
	      lineStart: lineStart,
	      lineEnd: lineEnd,
	      polygonStart: function() {
	        listener = bufferListener;
	        segments = [];
	        polygon = [];
	        clean = true;
	      },
	      polygonEnd: function() {
	        listener = listener_;
	        segments = d3f.merge(segments);
	        var clipStartInside = insidePolygon([x0, y1]),
	            inside = clean && clipStartInside,
	            visible = segments.length;
	        if (inside || visible) {
	          listener.polygonStart();
	          if (inside) {
	            listener.lineStart();
	            interpolate(null, null, 1, listener);
	            listener.lineEnd();
	          }
	          if (visible) {
	            d3f_geo_clipPolygon(segments, compare, clipStartInside, interpolate, listener);
	          }
	          listener.polygonEnd();
	        }
	        segments = polygon = ring = null;
	      }
	    };

	    function insidePolygon(p) {
	      var wn = 0, // the winding number counter
	          n = polygon.length,
	          y = p[1];

	      for (var i = 0; i < n; ++i) {
	        for (var j = 1, v = polygon[i], m = v.length, a = v[0], b; j < m; ++j) {
	          b = v[j];
	          if (a[1] <= y) {
	            if (b[1] >  y && d3f_cross2d(a, b, p) > 0) ++wn;
	          } else {
	            if (b[1] <= y && d3f_cross2d(a, b, p) < 0) --wn;
	          }
	          a = b;
	        }
	      }
	      return wn !== 0;
	    }

	    function interpolate(from, to, direction, listener) {
	      var a = 0, a1 = 0;
	      if (from == null ||
	          (a = corner(from, direction)) !== (a1 = corner(to, direction)) ||
	          comparePoints(from, to) < 0 ^ direction > 0) {
	        do {
	          listener.point(a === 0 || a === 3 ? x0 : x1, a > 1 ? y1 : y0);
	        } while ((a = (a + direction + 4) % 4) !== a1);
	      } else {
	        listener.point(to[0], to[1]);
	      }
	    }

	    function pointVisible(x, y) {
	      return x0 <= x && x <= x1 && y0 <= y && y <= y1;
	    }

	    function point(x, y) {
	      if (pointVisible(x, y)) listener.point(x, y);
	    }

	    var x__, y__, v__, // first point
	        x_, y_, v_, // previous point
	        first,
	        clean;

	    function lineStart() {
	      clip.point = linePoint;
	      if (polygon) polygon.push(ring = []);
	      first = true;
	      v_ = false;
	      x_ = y_ = NaN;
	    }

	    function lineEnd() {
	      // TODO rather than special-case polygons, simply handle them separately.
	      // Ideally, coincident intersection points should be jittered to avoid
	      // clipping issues.
	      if (segments) {
	        linePoint(x__, y__);
	        if (v__ && v_) bufferListener.rejoin();
	        segments.push(bufferListener.buffer());
	      }
	      clip.point = point;
	      if (v_) listener.lineEnd();
	    }

	    function linePoint(x, y) {
	      x = Math.max(-d3f_geo_clipExtentMAX, Math.min(d3f_geo_clipExtentMAX, x));
	      y = Math.max(-d3f_geo_clipExtentMAX, Math.min(d3f_geo_clipExtentMAX, y));
	      var v = pointVisible(x, y);
	      if (polygon) ring.push([x, y]);
	      if (first) {
	        x__ = x, y__ = y, v__ = v;
	        first = false;
	        if (v) {
	          listener.lineStart();
	          listener.point(x, y);
	        }
	      } else {
	        if (v && v_) listener.point(x, y);
	        else {
	          var l = {a: {x: x_, y: y_}, b: {x: x, y: y}};
	          if (clipLine(l)) {
	            if (!v_) {
	              listener.lineStart();
	              listener.point(l.a.x, l.a.y);
	            }
	            listener.point(l.b.x, l.b.y);
	            if (!v) listener.lineEnd();
	            clean = false;
	          } else if (v) {
	            listener.lineStart();
	            listener.point(x, y);
	            clean = false;
	          }
	        }
	      }
	      x_ = x, y_ = y, v_ = v;
	    }

	    return clip;
	  };

	  function corner(p, direction) {
	    return abs(p[0] - x0) < ε ? direction > 0 ? 0 : 3
	        : abs(p[0] - x1) < ε ? direction > 0 ? 2 : 1
	        : abs(p[1] - y0) < ε ? direction > 0 ? 1 : 0
	        : direction > 0 ? 3 : 2; // abs(p[1] - y1) < ε
	  }

	  function compare(a, b) {
	    return comparePoints(a.x, b.x);
	  }

	  function comparePoints(a, b) {
	    var ca = corner(a, 1),
	        cb = corner(b, 1);
	    return ca !== cb ? ca - cb
	        : ca === 0 ? b[1] - a[1]
	        : ca === 1 ? a[0] - b[0]
	        : ca === 2 ? a[1] - b[1]
	        : b[0] - a[0];
	  }
	}

	function d3f_geo_resample(project) {
	  var δ2 = .5, // precision, px²
	      cosMinDistance = Math.cos(30 * d3f_radians), // cos(minimum angular distance)
	      maxDepth = 16;

	  function resample(stream) {
	    return (maxDepth ? resampleRecursive : resampleNone)(stream);
	  }

	  function resampleNone(stream) {
	    return d3f_geo_transformPoint(stream, function(x, y) {
	      x = project(x, y);
	      stream.point(x[0], x[1]);
	    });
	  }

	  function resampleRecursive(stream) {
	    var λ00, φ00, x00, y00, a00, b00, c00, // first point
	        λ0, x0, y0, a0, b0, c0; // previous point

	    var resample = {
	      point: point,
	      lineStart: lineStart,
	      lineEnd: lineEnd,
	      polygonStart: function() { stream.polygonStart(); resample.lineStart = ringStart; },
	      polygonEnd: function() { stream.polygonEnd(); resample.lineStart = lineStart; }
	    };

	    function point(x, y) {
	      x = project(x, y);
	      stream.point(x[0], x[1]);
	    }

	    function lineStart() {
	      x0 = NaN;
	      resample.point = linePoint;
	      stream.lineStart();
	    }

	    function linePoint(λ, φ) {
	      var c = d3f_geo_cartesian([λ, φ]), p = project(λ, φ);
	      resampleLineTo(x0, y0, λ0, a0, b0, c0, x0 = p[0], y0 = p[1], λ0 = λ, a0 = c[0], b0 = c[1], c0 = c[2], maxDepth, stream);
	      stream.point(x0, y0);
	    }

	    function lineEnd() {
	      resample.point = point;
	      stream.lineEnd();
	    }

	    function ringStart() {
	      lineStart();
	      resample.point = ringPoint;
	      resample.lineEnd = ringEnd;
	    }

	    function ringPoint(λ, φ) {
	      linePoint(λ00 = λ, φ00 = φ), x00 = x0, y00 = y0, a00 = a0, b00 = b0, c00 = c0;
	      resample.point = linePoint;
	    }

	    function ringEnd() {
	      resampleLineTo(x0, y0, λ0, a0, b0, c0, x00, y00, λ00, a00, b00, c00, maxDepth, stream);
	      resample.lineEnd = lineEnd;
	      lineEnd();
	    }

	    return resample;
	  }

	  function resampleLineTo(x0, y0, λ0, a0, b0, c0, x1, y1, λ1, a1, b1, c1, depth, stream) {
	    var dx = x1 - x0,
	        dy = y1 - y0,
	        d2 = dx * dx + dy * dy;
	    if (d2 > 4 * δ2 && depth--) {
	      var a = a0 + a1,
	          b = b0 + b1,
	          c = c0 + c1,
	          m = Math.sqrt(a * a + b * b + c * c),
	          φ2 = Math.asin(c /= m),
	          λ2 = abs(abs(c) - 1) < ε || abs(λ0 - λ1) < ε ? (λ0 + λ1) / 2 : Math.atan2(b, a),
	          p = project(λ2, φ2),
	          x2 = p[0],
	          y2 = p[1],
	          dx2 = x2 - x0,
	          dy2 = y2 - y0,
	          dz = dy * dx2 - dx * dy2;
	      if (dz * dz / d2 > δ2 // perpendicular projected distance
	          || abs((dx * dx2 + dy * dy2) / d2 - .5) > .3 // midpoint close to an end
	          || a0 * a1 + b0 * b1 + c0 * c1 < cosMinDistance) { // angular distance
	        resampleLineTo(x0, y0, λ0, a0, b0, c0, x2, y2, λ2, a /= m, b /= m, c, depth, stream);
	        stream.point(x2, y2);
	        resampleLineTo(x2, y2, λ2, a, b, c, x1, y1, λ1, a1, b1, c1, depth, stream);
	      }
	    }
	  }

	  resample.precision = function(_) {
	    if (!arguments.length) return Math.sqrt(δ2);
	    maxDepth = (δ2 = _ * _) > 0 && 16;
	    return resample;
	  };

	  return resample;
	}
	var d3f_arraySlice = [].slice,
	    d3f_array = function(list) { return d3f_arraySlice.call(list); }; // conversion for NodeLists

	d3f.geo.transform = function(methods) {
	  return {
	    stream: function(stream) {
	      var transform = new d3f_geo_transform(stream);
	      for (var k in methods) transform[k] = methods[k];
	      return transform;
	    }
	  };
	};

	function d3f_geo_transform(stream) {
	  this.stream = stream;
	}

	d3f_geo_transform.prototype = {
	  point: function(x, y) { this.stream.point(x, y); },
	  sphere: function() { this.stream.sphere(); },
	  lineStart: function() { this.stream.lineStart(); },
	  lineEnd: function() { this.stream.lineEnd(); },
	  polygonStart: function() { this.stream.polygonStart(); },
	  polygonEnd: function() { this.stream.polygonEnd(); }
	};

	function d3f_geo_transformPoint(stream, point) {
	  return {
	    point: point,
	    sphere: function() { stream.sphere(); },
	    lineStart: function() { stream.lineStart(); },
	    lineEnd: function() { stream.lineEnd(); },
	    polygonStart: function() { stream.polygonStart(); },
	    polygonEnd: function() { stream.polygonEnd(); },
	  };
	}

	d3f.geo.projection = d3f_geo_projection;
	d3f.geo.projectionMutator = d3f_geo_projectionMutator;

	function d3f_geo_projection(project) {
	  return d3f_geo_projectionMutator(function() { return project; })();
	}

	function d3f_geo_projectionMutator(projectAt) {
	  var project,
	      rotate,
	      projectRotate,
	      projectResample = d3f_geo_resample(function(x, y) { x = project(x, y); return [x[0] * k + δx, δy - x[1] * k]; }),
	      k = 150, // scale
	      x = 480, y = 250, // translate
	      λ = 0, φ = 0, // center
	      δλ = 0, δφ = 0, δγ = 0, // rotate
	      δx, δy, // center
	      preclip = d3f_geo_clipAntimeridian,
	      postclip = d3f_identity,
	      clipAngle = null,
	      clipExtent = null,
	      stream;

	  function projection(point) {
	    point = projectRotate(point[0] * d3f_radians, point[1] * d3f_radians);
	    return [point[0] * k + δx, δy - point[1] * k];
	  }

	  function invert(point) {
	    point = projectRotate.invert((point[0] - δx) / k, (δy - point[1]) / k);
	    return point && [point[0] * d3f_degrees, point[1] * d3f_degrees];
	  }

	  projection.stream = function(output) {
	    if (stream) stream.valid = false;
	    stream = d3f_geo_projectionRadians(preclip(rotate, projectResample(postclip(output))));
	    stream.valid = true; // allow caching by d3f.geo.path
	    return stream;
	  };

	  projection.clipAngle = function(_) {
	    if (!arguments.length) return clipAngle;
	    preclip = _ == null ? (clipAngle = _, d3f_geo_clipAntimeridian) : d3f_geo_clipCircle((clipAngle = +_) * d3f_radians);
	    return invalidate();
	  };

	  projection.clipExtent = function(_) {
	    if (!arguments.length) return clipExtent;
	    clipExtent = _;
	    postclip = _ ? d3f_geo_clipExtent(_[0][0], _[0][1], _[1][0], _[1][1]) : d3f_identity;
	    return invalidate();
	  };

	  projection.scale = function(_) {
	    if (!arguments.length) return k;
	    k = +_;
	    return reset();
	  };

	  projection.translate = function(_) {
	    if (!arguments.length) return [x, y];
	    x = +_[0];
	    y = +_[1];
	    return reset();
	  };

	  projection.center = function(_) {
	    if (!arguments.length) return [λ * d3f_degrees, φ * d3f_degrees];
	    λ = _[0] % 360 * d3f_radians;
	    φ = _[1] % 360 * d3f_radians;
	    return reset();
	  };

	  projection.rotate = function(_) {
	    if (!arguments.length) return [δλ * d3f_degrees, δφ * d3f_degrees, δγ * d3f_degrees];
	    δλ = _[0] % 360 * d3f_radians;
	    δφ = _[1] % 360 * d3f_radians;
	    δγ = _.length > 2 ? _[2] % 360 * d3f_radians : 0;
	    return reset();
	  };

	  d3f.rebind(projection, projectResample, "precision");

	  function reset() {
	    projectRotate = d3f_geo_compose(rotate = d3f_geo_rotation(δλ, δφ, δγ), project);
	    var center = project(λ, φ);
	    δx = x - center[0] * k;
	    δy = y + center[1] * k;
	    return invalidate();
	  }

	  function invalidate() {
	    if (stream) stream.valid = false, stream = null;
	    return projection;
	  }

	  return function() {
	    project = projectAt.apply(this, arguments);
	    projection.invert = project.invert && invert;
	    return reset();
	  };
	}

	function d3f_geo_projectionRadians(stream) {
	  return d3f_geo_transformPoint(stream, function(x, y) {
	    stream.point(x * d3f_radians, y * d3f_radians);
	  });
	}

	function d3f_geo_conic(projectAt) {
	  var φ0 = 0,
	      φ1 = π / 3,
	      m = d3f_geo_projectionMutator(projectAt),
	      p = m(φ0, φ1);

	  p.parallels = function(_) {
	    if (!arguments.length) return [φ0 / π * 180, φ1 / π * 180];
	    return m(φ0 = _[0] * π / 180, φ1 = _[1] * π / 180);
	  };

	  return p;
	}

	function d3f_geo_conicEqualArea(φ0, φ1) {
	  var sinφ0 = Math.sin(φ0),
	      n = (sinφ0 + Math.sin(φ1)) / 2,
	      C = 1 + sinφ0 * (2 * n - sinφ0),
	      ρ0 = Math.sqrt(C) / n;

	  function forward(λ, φ) {
	    var ρ = Math.sqrt(C - 2 * n * Math.sin(φ)) / n;
	    return [
	      ρ * Math.sin(λ *= n),
	      ρ0 - ρ * Math.cos(λ)
	    ];
	  }

	  forward.invert = function(x, y) {
	    var ρ0_y = ρ0 - y;
	    return [
	      Math.atan2(x, ρ0_y) / n,
	      d3f_asin((C - (x * x + ρ0_y * ρ0_y) * n * n) / (2 * n))
	    ];
	  };

	  return forward;
	}

	(d3f.geo.conicEqualArea = function() {
	  return d3f_geo_conic(d3f_geo_conicEqualArea);
	}).raw = d3f_geo_conicEqualArea;

	// ESRI:102003
	d3f.geo.albers = function() {
	  return d3f.geo.conicEqualArea()
	      .rotate([96, 0])
	      .center([-.6, 38.7])
	      .parallels([29.5, 45.5])
	      .scale(1070);
	};

	// A composite projection for the United States, configured by default for
	// 960×500. Also works quite well at 960×600 with scale 1285. The set of
	// standard parallels for each region comes from USGS, which is published here:
	// http://egsc.usgs.gov/isb/pubs/MapProjections/projections.html#albers
	d3f.geo.albersUsa = function() {
	  var lower48 = d3f.geo.albers();

	  // EPSG:3338
	  var alaska = d3f.geo.conicEqualArea()
	      .rotate([154, 0])
	      .center([-2, 58.5])
	      .parallels([55, 65]);

	  // ESRI:102007
	  var hawaii = d3f.geo.conicEqualArea()
	      .rotate([157, 0])
	      .center([-3, 19.9])
	      .parallels([8, 18]);

	  var point,
	      pointStream = {point: function(x, y) { point = [x, y]; }},
	      lower48Point,
	      alaskaPoint,
	      hawaiiPoint;

	  function albersUsa(coordinates) {
	    var x = coordinates[0], y = coordinates[1];
	    point = null;
	    (lower48Point(x, y), point)
	        || (alaskaPoint(x, y), point)
	        || hawaiiPoint(x, y);
	    return point;
	  }

	  albersUsa.invert = function(coordinates) {
	    var k = lower48.scale(),
	        t = lower48.translate(),
	        x = (coordinates[0] - t[0]) / k,
	        y = (coordinates[1] - t[1]) / k;
	    return (y >= .120 && y < .234 && x >= -.425 && x < -.214 ? alaska
	        : y >= .166 && y < .234 && x >= -.214 && x < -.115 ? hawaii
	        : lower48).invert(coordinates);
	  };

	  // A naïve multi-projection stream.
	  // The projections must have mutually exclusive clip regions on the sphere,
	  // as this will avoid emitting interleaving lines and polygons.
	  albersUsa.stream = function(stream) {
	    var lower48Stream = lower48.stream(stream),
	        alaskaStream = alaska.stream(stream),
	        hawaiiStream = hawaii.stream(stream);
	    return {
	      point: function(x, y) {
	        lower48Stream.point(x, y);
	        alaskaStream.point(x, y);
	        hawaiiStream.point(x, y);
	      },
	      sphere: function() {
	        lower48Stream.sphere();
	        alaskaStream.sphere();
	        hawaiiStream.sphere();
	      },
	      lineStart: function() {
	        lower48Stream.lineStart();
	        alaskaStream.lineStart();
	        hawaiiStream.lineStart();
	      },
	      lineEnd: function() {
	        lower48Stream.lineEnd();
	        alaskaStream.lineEnd();
	        hawaiiStream.lineEnd();
	      },
	      polygonStart: function() {
	        lower48Stream.polygonStart();
	        alaskaStream.polygonStart();
	        hawaiiStream.polygonStart();
	      },
	      polygonEnd: function() {
	        lower48Stream.polygonEnd();
	        alaskaStream.polygonEnd();
	        hawaiiStream.polygonEnd();
	      }
	    };
	  };

	  albersUsa.precision = function(_) {
	    if (!arguments.length) return lower48.precision();
	    lower48.precision(_);
	    alaska.precision(_);
	    hawaii.precision(_);
	    return albersUsa;
	  };

	  albersUsa.scale = function(_) {
	    if (!arguments.length) return lower48.scale();
	    lower48.scale(_);
	    alaska.scale(_ * .35);
	    hawaii.scale(_);
	    return albersUsa.translate(lower48.translate());
	  };

	  albersUsa.translate = function(_) {
	    if (!arguments.length) return lower48.translate();
	    var k = lower48.scale(), x = +_[0], y = +_[1];

	    lower48Point = lower48
	        .translate(_)
	        .clipExtent([[x - .455 * k, y - .238 * k], [x + .455 * k, y + .238 * k]])
	        .stream(pointStream).point;

	    alaskaPoint = alaska
	        .translate([x - .307 * k, y + .201 * k])
	        .clipExtent([[x - .425 * k + ε, y + .120 * k + ε], [x - .214 * k - ε, y + .234 * k - ε]])
	        .stream(pointStream).point;

	    hawaiiPoint = hawaii
	        .translate([x - .205 * k, y + .212 * k])
	        .clipExtent([[x - .214 * k + ε, y + .166 * k + ε], [x - .115 * k - ε, y + .234 * k - ε]])
	        .stream(pointStream).point;

	    return albersUsa;
	  };

	  return albersUsa.scale(1070);
	};

	d3f.geo.bounds = (function() {
	  var λ0, φ0, λ1, φ1, // bounds
	      λ_, // previous λ-coordinate
	      λ__, φ__, // first point
	      p0, // previous 3D point
	      dλSum,
	      ranges,
	      range;

	  var bound = {
	    point: point,
	    lineStart: lineStart,
	    lineEnd: lineEnd,

	    polygonStart: function() {
	      bound.point = ringPoint;
	      bound.lineStart = ringStart;
	      bound.lineEnd = ringEnd;
	      dλSum = 0;
	      d3f_geo_area.polygonStart();
	    },
	    polygonEnd: function() {
	      d3f_geo_area.polygonEnd();
	      bound.point = point;
	      bound.lineStart = lineStart;
	      bound.lineEnd = lineEnd;
	      if (d3f_geo_areaRingSum < 0) λ0 = -(λ1 = 180), φ0 = -(φ1 = 90);
	      else if (dλSum > ε) φ1 = 90;
	      else if (dλSum < -ε) φ0 = -90;
	      range[0] = λ0, range[1] = λ1;
	    }
	  };

	  function point(λ, φ) {
	    ranges.push(range = [λ0 = λ, λ1 = λ]);
	    if (φ < φ0) φ0 = φ;
	    if (φ > φ1) φ1 = φ;
	  }

	  function linePoint(λ, φ) {
	    var p = d3f_geo_cartesian([λ * d3f_radians, φ * d3f_radians]);
	    if (p0) {
	      var normal = d3f_geo_cartesianCross(p0, p),
	          equatorial = [normal[1], -normal[0], 0],
	          inflection = d3f_geo_cartesianCross(equatorial, normal);
	      d3f_geo_cartesianNormalize(inflection);
	      inflection = d3f_geo_spherical(inflection);
	      var dλ = λ - λ_,
	          s = dλ > 0 ? 1 : -1,
	          λi = inflection[0] * d3f_degrees * s,
	          antimeridian = abs(dλ) > 180;
	      if (antimeridian ^ (s * λ_ < λi && λi < s * λ)) {
	        var φi = inflection[1] * d3f_degrees;
	        if (φi > φ1) φ1 = φi;
	      } else if (λi = (λi + 360) % 360 - 180, antimeridian ^ (s * λ_ < λi && λi < s * λ)) {
	        var φi = -inflection[1] * d3f_degrees;
	        if (φi < φ0) φ0 = φi;
	      } else {
	        if (φ < φ0) φ0 = φ;
	        if (φ > φ1) φ1 = φ;
	      }
	      if (antimeridian) {
	        if (λ < λ_) {
	          if (angle(λ0, λ) > angle(λ0, λ1)) λ1 = λ;
	        } else {
	          if (angle(λ, λ1) > angle(λ0, λ1)) λ0 = λ;
	        }
	      } else {
	        if (λ1 >= λ0) {
	          if (λ < λ0) λ0 = λ;
	          if (λ > λ1) λ1 = λ;
	        } else {
	          if (λ > λ_) {
	            if (angle(λ0, λ) > angle(λ0, λ1)) λ1 = λ;
	          } else {
	            if (angle(λ, λ1) > angle(λ0, λ1)) λ0 = λ;
	          }
	        }
	      }
	    } else {
	      point(λ, φ);
	    }
	    p0 = p, λ_ = λ;
	  }

	  function lineStart() { bound.point = linePoint; }
	  function lineEnd() {
	    range[0] = λ0, range[1] = λ1;
	    bound.point = point;
	    p0 = null;
	  }

	  function ringPoint(λ, φ) {
	    if (p0) {
	      var dλ = λ - λ_;
	      dλSum += abs(dλ) > 180 ? dλ + (dλ > 0 ? 360 : -360) : dλ;
	    } else λ__ = λ, φ__ = φ;
	    d3f_geo_area.point(λ, φ);
	    linePoint(λ, φ);
	  }

	  function ringStart() {
	    d3f_geo_area.lineStart();
	  }

	  function ringEnd() {
	    ringPoint(λ__, φ__);
	    d3f_geo_area.lineEnd();
	    if (abs(dλSum) > ε) λ0 = -(λ1 = 180);
	    range[0] = λ0, range[1] = λ1;
	    p0 = null;
	  }

	  // Finds the left-right distance between two longitudes.
	  // This is almost the same as (λ1 - λ0 + 360°) % 360°, except that we want
	  // the distance between ±180° to be 360°.
	  function angle(λ0, λ1) { return (λ1 -= λ0) < 0 ? λ1 + 360 : λ1; }

	  function compareRanges(a, b) { return a[0] - b[0]; }

	  function withinRange(x, range) {
	    return range[0] <= range[1] ? range[0] <= x && x <= range[1] : x < range[0] || range[1] < x;
	  }

	  return function(feature) {
	    φ1 = λ1 = -(λ0 = φ0 = Infinity);
	    ranges = [];

	    d3f.geo.stream(feature, bound);

	    var n = ranges.length;
	    if (n) {
	      // First, sort ranges by their minimum longitudes.
	      ranges.sort(compareRanges);

	      // Then, merge any ranges that overlap.
	      for (var i = 1, a = ranges[0], b, merged = [a]; i < n; ++i) {
	        b = ranges[i];
	        if (withinRange(b[0], a) || withinRange(b[1], a)) {
	          if (angle(a[0], b[1]) > angle(a[0], a[1])) a[1] = b[1];
	          if (angle(b[0], a[1]) > angle(a[0], a[1])) a[0] = b[0];
	        } else {
	          merged.push(a = b);
	        }
	      }

	      // Finally, find the largest gap between the merged ranges.
	      // The final bounding box will be the inverse of this gap.
	      var best = -Infinity, dλ;
	      for (var n = merged.length - 1, i = 0, a = merged[n], b; i <= n; a = b, ++i) {
	        b = merged[i];
	        if ((dλ = angle(a[1], b[0])) > best) best = dλ, λ0 = b[0], λ1 = a[1];
	      }
	    }
	    ranges = range = null;

	    return λ0 === Infinity || φ0 === Infinity
	        ? [[NaN, NaN], [NaN, NaN]]
	        : [[λ0, φ0], [λ1, φ1]];
	  };
	})();

	d3f.geo.centroid = function(object) {
	  d3f_geo_centroidW0 = d3f_geo_centroidW1 =
	  d3f_geo_centroidX0 = d3f_geo_centroidY0 = d3f_geo_centroidZ0 =
	  d3f_geo_centroidX1 = d3f_geo_centroidY1 = d3f_geo_centroidZ1 =
	  d3f_geo_centroidX2 = d3f_geo_centroidY2 = d3f_geo_centroidZ2 = 0;
	  d3f.geo.stream(object, d3f_geo_centroid);

	  var x = d3f_geo_centroidX2,
	      y = d3f_geo_centroidY2,
	      z = d3f_geo_centroidZ2,
	      m = x * x + y * y + z * z;

	  // If the area-weighted centroid is undefined, fall back to length-weighted centroid.
	  if (m < ε2) {
	    x = d3f_geo_centroidX1, y = d3f_geo_centroidY1, z = d3f_geo_centroidZ1;
	    // If the feature has zero length, fall back to arithmetic mean of point vectors.
	    if (d3f_geo_centroidW1 < ε) x = d3f_geo_centroidX0, y = d3f_geo_centroidY0, z = d3f_geo_centroidZ0;
	    m = x * x + y * y + z * z;
	    // If the feature still has an undefined centroid, then return.
	    if (m < ε2) return [NaN, NaN];
	  }

	  return [Math.atan2(y, x) * d3f_degrees, d3f_asin(z / Math.sqrt(m)) * d3f_degrees];
	};

	var d3f_geo_centroidW0,
	    d3f_geo_centroidW1,
	    d3f_geo_centroidX0,
	    d3f_geo_centroidY0,
	    d3f_geo_centroidZ0,
	    d3f_geo_centroidX1,
	    d3f_geo_centroidY1,
	    d3f_geo_centroidZ1,
	    d3f_geo_centroidX2,
	    d3f_geo_centroidY2,
	    d3f_geo_centroidZ2;

	var d3f_geo_centroid = {
	  sphere: d3f_noop,
	  point: d3f_geo_centroidPoint,
	  lineStart: d3f_geo_centroidLineStart,
	  lineEnd: d3f_geo_centroidLineEnd,
	  polygonStart: function() {
	    d3f_geo_centroid.lineStart = d3f_geo_centroidRingStart;
	  },
	  polygonEnd: function() {
	    d3f_geo_centroid.lineStart = d3f_geo_centroidLineStart;
	  }
	};

	// Arithmetic mean of Cartesian vectors.
	function d3f_geo_centroidPoint(λ, φ) {
	  λ *= d3f_radians;
	  var cosφ = Math.cos(φ *= d3f_radians);
	  d3f_geo_centroidPointXYZ(cosφ * Math.cos(λ), cosφ * Math.sin(λ), Math.sin(φ));
	}

	function d3f_geo_centroidPointXYZ(x, y, z) {
	  ++d3f_geo_centroidW0;
	  d3f_geo_centroidX0 += (x - d3f_geo_centroidX0) / d3f_geo_centroidW0;
	  d3f_geo_centroidY0 += (y - d3f_geo_centroidY0) / d3f_geo_centroidW0;
	  d3f_geo_centroidZ0 += (z - d3f_geo_centroidZ0) / d3f_geo_centroidW0;
	}

	function d3f_geo_centroidLineStart() {
	  var x0, y0, z0; // previous point

	  d3f_geo_centroid.point = function(λ, φ) {
	    λ *= d3f_radians;
	    var cosφ = Math.cos(φ *= d3f_radians);
	    x0 = cosφ * Math.cos(λ);
	    y0 = cosφ * Math.sin(λ);
	    z0 = Math.sin(φ);
	    d3f_geo_centroid.point = nextPoint;
	    d3f_geo_centroidPointXYZ(x0, y0, z0);
	  };

	  function nextPoint(λ, φ) {
	    λ *= d3f_radians;
	    var cosφ = Math.cos(φ *= d3f_radians),
	        x = cosφ * Math.cos(λ),
	        y = cosφ * Math.sin(λ),
	        z = Math.sin(φ),
	        w = Math.atan2(
	          Math.sqrt((w = y0 * z - z0 * y) * w + (w = z0 * x - x0 * z) * w + (w = x0 * y - y0 * x) * w),
	          x0 * x + y0 * y + z0 * z);
	    d3f_geo_centroidW1 += w;
	    d3f_geo_centroidX1 += w * (x0 + (x0 = x));
	    d3f_geo_centroidY1 += w * (y0 + (y0 = y));
	    d3f_geo_centroidZ1 += w * (z0 + (z0 = z));
	    d3f_geo_centroidPointXYZ(x0, y0, z0);
	  }
	}

	function d3f_geo_centroidLineEnd() {
	  d3f_geo_centroid.point = d3f_geo_centroidPoint;
	}

	// See J. E. Brock, The Inertia Tensor for a Spherical Triangle,
	// J. Applied Mechanics 42, 239 (1975).
	function d3f_geo_centroidRingStart() {
	  var λ00, φ00, // first point
	      x0, y0, z0; // previous point

	  d3f_geo_centroid.point = function(λ, φ) {
	    λ00 = λ, φ00 = φ;
	    d3f_geo_centroid.point = nextPoint;
	    λ *= d3f_radians;
	    var cosφ = Math.cos(φ *= d3f_radians);
	    x0 = cosφ * Math.cos(λ);
	    y0 = cosφ * Math.sin(λ);
	    z0 = Math.sin(φ);
	    d3f_geo_centroidPointXYZ(x0, y0, z0);
	  };

	  d3f_geo_centroid.lineEnd = function() {
	    nextPoint(λ00, φ00);
	    d3f_geo_centroid.lineEnd = d3f_geo_centroidLineEnd;
	    d3f_geo_centroid.point = d3f_geo_centroidPoint;
	  };

	  function nextPoint(λ, φ) {
	    λ *= d3f_radians;
	    var cosφ = Math.cos(φ *= d3f_radians),
	        x = cosφ * Math.cos(λ),
	        y = cosφ * Math.sin(λ),
	        z = Math.sin(φ),
	        cx = y0 * z - z0 * y,
	        cy = z0 * x - x0 * z,
	        cz = x0 * y - y0 * x,
	        m = Math.sqrt(cx * cx + cy * cy + cz * cz),
	        u = x0 * x + y0 * y + z0 * z,
	        v = m && -d3f_acos(u) / m, // area weight
	        w = Math.atan2(m, u); // line weight
	    d3f_geo_centroidX2 += v * cx;
	    d3f_geo_centroidY2 += v * cy;
	    d3f_geo_centroidZ2 += v * cz;
	    d3f_geo_centroidW1 += w;
	    d3f_geo_centroidX1 += w * (x0 + (x0 = x));
	    d3f_geo_centroidY1 += w * (y0 + (y0 = y));
	    d3f_geo_centroidZ1 += w * (z0 + (z0 = z));
	    d3f_geo_centroidPointXYZ(x0, y0, z0);
	  }
	}

	// TODO Unify this code with d3f.geom.polygon area?

	var d3f_geo_pathAreaSum, d3f_geo_pathAreaPolygon, d3f_geo_pathArea = {
	  point: d3f_noop,
	  lineStart: d3f_noop,
	  lineEnd: d3f_noop,

	  // Only count area for polygon rings.
	  polygonStart: function() {
	    d3f_geo_pathAreaPolygon = 0;
	    d3f_geo_pathArea.lineStart = d3f_geo_pathAreaRingStart;
	  },
	  polygonEnd: function() {
	    d3f_geo_pathArea.lineStart = d3f_geo_pathArea.lineEnd = d3f_geo_pathArea.point = d3f_noop;
	    d3f_geo_pathAreaSum += abs(d3f_geo_pathAreaPolygon / 2);
	  }
	};

	function d3f_geo_pathAreaRingStart() {
	  var x00, y00, x0, y0;

	  // For the first point, …
	  d3f_geo_pathArea.point = function(x, y) {
	    d3f_geo_pathArea.point = nextPoint;
	    x00 = x0 = x, y00 = y0 = y;
	  };

	  // For subsequent points, …
	  function nextPoint(x, y) {
	    d3f_geo_pathAreaPolygon += y0 * x - x0 * y;
	    x0 = x, y0 = y;
	  }

	  // For the last point, return to the start.
	  d3f_geo_pathArea.lineEnd = function() {
	    nextPoint(x00, y00);
	  };
	}

	var d3f_geo_pathBoundsX0,
	    d3f_geo_pathBoundsY0,
	    d3f_geo_pathBoundsX1,
	    d3f_geo_pathBoundsY1;

	var d3f_geo_pathBounds = {
	  point: d3f_geo_pathBoundsPoint,
	  lineStart: d3f_noop,
	  lineEnd: d3f_noop,
	  polygonStart: d3f_noop,
	  polygonEnd: d3f_noop
	};

	function d3f_geo_pathBoundsPoint(x, y) {
	  if (x < d3f_geo_pathBoundsX0) d3f_geo_pathBoundsX0 = x;
	  if (x > d3f_geo_pathBoundsX1) d3f_geo_pathBoundsX1 = x;
	  if (y < d3f_geo_pathBoundsY0) d3f_geo_pathBoundsY0 = y;
	  if (y > d3f_geo_pathBoundsY1) d3f_geo_pathBoundsY1 = y;
	}
	function d3f_geo_pathBuffer() {
	  var pointCircle = d3f_geo_pathBufferCircle(4.5),
	      buffer = [];

	  var stream = {
	    point: point,

	    // While inside a line, override point to moveTo then lineTo.
	    lineStart: function() { stream.point = pointLineStart; },
	    lineEnd: lineEnd,

	    // While inside a polygon, override lineEnd to closePath.
	    polygonStart: function() { stream.lineEnd = lineEndPolygon; },
	    polygonEnd: function() { stream.lineEnd = lineEnd; stream.point = point; },

	    pointRadius: function(_) {
	      pointCircle = d3f_geo_pathBufferCircle(_);
	      return stream;
	    },

	    result: function() {
	      if (buffer.length) {
	        var result = buffer.join("");
	        buffer = [];
	        return result;
	      }
	    }
	  };

	  function point(x, y) {
	    buffer.push("M", x, ",", y, pointCircle);
	  }

	  function pointLineStart(x, y) {
	    buffer.push("M", x, ",", y);
	    stream.point = pointLine;
	  }

	  function pointLine(x, y) {
	    buffer.push("L", x, ",", y);
	  }

	  function lineEnd() {
	    stream.point = point;
	  }

	  function lineEndPolygon() {
	    buffer.push("Z");
	  }

	  return stream;
	}

	function d3f_geo_pathBufferCircle(radius) {
	  return "m0," + radius
	      + "a" + radius + "," + radius + " 0 1,1 0," + -2 * radius
	      + "a" + radius + "," + radius + " 0 1,1 0," + 2 * radius
	      + "z";
	}

	// TODO Unify this code with d3f.geom.polygon centroid?
	// TODO Enforce positive area for exterior, negative area for interior?

	var d3f_geo_pathCentroid = {
	  point: d3f_geo_pathCentroidPoint,

	  // For lines, weight by length.
	  lineStart: d3f_geo_pathCentroidLineStart,
	  lineEnd: d3f_geo_pathCentroidLineEnd,

	  // For polygons, weight by area.
	  polygonStart: function() {
	    d3f_geo_pathCentroid.lineStart = d3f_geo_pathCentroidRingStart;
	  },
	  polygonEnd: function() {
	    d3f_geo_pathCentroid.point = d3f_geo_pathCentroidPoint;
	    d3f_geo_pathCentroid.lineStart = d3f_geo_pathCentroidLineStart;
	    d3f_geo_pathCentroid.lineEnd = d3f_geo_pathCentroidLineEnd;
	  }
	};

	function d3f_geo_pathCentroidPoint(x, y) {
	  d3f_geo_centroidX0 += x;
	  d3f_geo_centroidY0 += y;
	  ++d3f_geo_centroidZ0;
	}

	function d3f_geo_pathCentroidLineStart() {
	  var x0, y0;

	  d3f_geo_pathCentroid.point = function(x, y) {
	    d3f_geo_pathCentroid.point = nextPoint;
	    d3f_geo_pathCentroidPoint(x0 = x, y0 = y);
	  };

	  function nextPoint(x, y) {
	    var dx = x - x0, dy = y - y0, z = Math.sqrt(dx * dx + dy * dy);
	    d3f_geo_centroidX1 += z * (x0 + x) / 2;
	    d3f_geo_centroidY1 += z * (y0 + y) / 2;
	    d3f_geo_centroidZ1 += z;
	    d3f_geo_pathCentroidPoint(x0 = x, y0 = y);
	  }
	}

	function d3f_geo_pathCentroidLineEnd() {
	  d3f_geo_pathCentroid.point = d3f_geo_pathCentroidPoint;
	}

	function d3f_geo_pathCentroidRingStart() {
	  var x00, y00, x0, y0;

	  // For the first point, …
	  d3f_geo_pathCentroid.point = function(x, y) {
	    d3f_geo_pathCentroid.point = nextPoint;
	    d3f_geo_pathCentroidPoint(x00 = x0 = x, y00 = y0 = y);
	  };

	  // For subsequent points, …
	  function nextPoint(x, y) {
	    var dx = x - x0, dy = y - y0, z = Math.sqrt(dx * dx + dy * dy);
	    d3f_geo_centroidX1 += z * (x0 + x) / 2;
	    d3f_geo_centroidY1 += z * (y0 + y) / 2;
	    d3f_geo_centroidZ1 += z;

	    z = y0 * x - x0 * y;
	    d3f_geo_centroidX2 += z * (x0 + x);
	    d3f_geo_centroidY2 += z * (y0 + y);
	    d3f_geo_centroidZ2 += z * 3;
	    d3f_geo_pathCentroidPoint(x0 = x, y0 = y);
	  }

	  // For the last point, return to the start.
	  d3f_geo_pathCentroid.lineEnd = function() {
	    nextPoint(x00, y00);
	  };
	}

	function d3f_geo_pathContext(context) {
	  var pointRadius = 4.5;

	  var stream = {
	    point: point,

	    // While inside a line, override point to moveTo then lineTo.
	    lineStart: function() { stream.point = pointLineStart; },
	    lineEnd: lineEnd,

	    // While inside a polygon, override lineEnd to closePath.
	    polygonStart: function() { stream.lineEnd = lineEndPolygon; },
	    polygonEnd: function() { stream.lineEnd = lineEnd; stream.point = point; },

	    pointRadius: function(_) {
	      pointRadius = _;
	      return stream;
	    },

	    result: d3f_noop
	  };

	  function point(x, y) {
	    context.moveTo(x + pointRadius, y);
	    context.arc(x, y, pointRadius, 0, τ);
	  }

	  function pointLineStart(x, y) {
	    context.moveTo(x, y);
	    stream.point = pointLine;
	  }

	  function pointLine(x, y) {
	    context.lineTo(x, y);
	  }

	  function lineEnd() {
	    stream.point = point;
	  }

	  function lineEndPolygon() {
	    context.closePath();
	  }

	  return stream;
	}

	d3f.geo.path = function() {
	  var pointRadius = 4.5,
	      projection,
	      context,
	      projectStream,
	      contextStream,
	      cacheStream;

	  function path(object) {
	    if (object) {
	      if (typeof pointRadius === "function") contextStream.pointRadius(+pointRadius.apply(this, arguments));
	      if (!cacheStream || !cacheStream.valid) cacheStream = projectStream(contextStream);
	      d3f.geo.stream(object, cacheStream);
	    }
	    return contextStream.result();
	  }

	  path.area = function(object) {
	    d3f_geo_pathAreaSum = 0;
	    d3f.geo.stream(object, projectStream(d3f_geo_pathArea));
	    return d3f_geo_pathAreaSum;
	  };

	  path.centroid = function(object) {
	    d3f_geo_centroidX0 = d3f_geo_centroidY0 = d3f_geo_centroidZ0 =
	    d3f_geo_centroidX1 = d3f_geo_centroidY1 = d3f_geo_centroidZ1 =
	    d3f_geo_centroidX2 = d3f_geo_centroidY2 = d3f_geo_centroidZ2 = 0;
	    d3f.geo.stream(object, projectStream(d3f_geo_pathCentroid));
	    return d3f_geo_centroidZ2 ? [d3f_geo_centroidX2 / d3f_geo_centroidZ2, d3f_geo_centroidY2 / d3f_geo_centroidZ2]
	        : d3f_geo_centroidZ1 ? [d3f_geo_centroidX1 / d3f_geo_centroidZ1, d3f_geo_centroidY1 / d3f_geo_centroidZ1]
	        : d3f_geo_centroidZ0 ? [d3f_geo_centroidX0 / d3f_geo_centroidZ0, d3f_geo_centroidY0 / d3f_geo_centroidZ0]
	        : [NaN, NaN];
	  };

	  path.bounds = function(object) {
	    d3f_geo_pathBoundsX1 = d3f_geo_pathBoundsY1 = -(d3f_geo_pathBoundsX0 = d3f_geo_pathBoundsY0 = Infinity);
	    d3f.geo.stream(object, projectStream(d3f_geo_pathBounds));
	    return [[d3f_geo_pathBoundsX0, d3f_geo_pathBoundsY0], [d3f_geo_pathBoundsX1, d3f_geo_pathBoundsY1]];
	  };

	  path.projection = function(_) {
	    if (!arguments.length) return projection;
	    projectStream = (projection = _) ? _.stream || d3f_geo_pathProjectStream(_) : d3f_identity;
	    return reset();
	  };

	  path.context = function(_) {
	    if (!arguments.length) return context;
	    contextStream = (context = _) == null ? new d3f_geo_pathBuffer : new d3f_geo_pathContext(_);
	    if (typeof pointRadius !== "function") contextStream.pointRadius(pointRadius);
	    return reset();
	  };

	  path.pointRadius = function(_) {
	    if (!arguments.length) return pointRadius;
	    pointRadius = typeof _ === "function" ? _ : (contextStream.pointRadius(+_), +_);
	    return path;
	  };

	  function reset() {
	    cacheStream = null;
	    return path;
	  }

	  return path.projection(d3f.geo.albersUsa()).context(null);
	};

	function d3f_geo_pathProjectStream(project) {
	  var resample = d3f_geo_resample(function(x, y) { return project([x * d3f_degrees, y * d3f_degrees]); });
	  return function(stream) { return d3f_geo_projectionRadians(resample(stream)); };
	}

	function d3f_geo_mercator(λ, φ) {
	  return [λ, Math.log(Math.tan(π / 4 + φ / 2))];
	}

	d3f_geo_mercator.invert = function(x, y) {
	  return [x, 2 * Math.atan(Math.exp(y)) - halfπ];
	};

	function d3f_geo_mercatorProjection(project) {
	  var m = d3f_geo_projection(project),
	      scale = m.scale,
	      translate = m.translate,
	      clipExtent = m.clipExtent,
	      clipAuto;

	  m.scale = function() {
	    var v = scale.apply(m, arguments);
	    return v === m ? (clipAuto ? m.clipExtent(null) : m) : v;
	  };

	  m.translate = function() {
	    var v = translate.apply(m, arguments);
	    return v === m ? (clipAuto ? m.clipExtent(null) : m) : v;
	  };

	  m.clipExtent = function(_) {
	    var v = clipExtent.apply(m, arguments);
	    if (v === m) {
	      if (clipAuto = _ == null) {
	        var k = π * scale(), t = translate();
	        clipExtent([[t[0] - k, t[1] - k], [t[0] + k, t[1] + k]]);
	      }
	    } else if (clipAuto) {
	      v = null;
	    }
	    return v;
	  };

	  return m.clipExtent(null);
	}

	(d3f.geo.mercator = function() {
	  return d3f_geo_mercatorProjection(d3f_geo_mercator);
	}).raw = d3f_geo_mercator;
	  if (true) !(__WEBPACK_AMD_DEFINE_FACTORY__ = (d3f), __WEBPACK_AMD_DEFINE_RESULT__ = (typeof __WEBPACK_AMD_DEFINE_FACTORY__ === 'function' ? (__WEBPACK_AMD_DEFINE_FACTORY__.call(exports, __webpack_require__, exports, module)) : __WEBPACK_AMD_DEFINE_FACTORY__), __WEBPACK_AMD_DEFINE_RESULT__ !== undefined && (module.exports = __WEBPACK_AMD_DEFINE_RESULT__));
	  else if (typeof module === "object" && module.exports) module.exports = d3f;
	  this.d3f = d3f;
	}();


/***/ },
/* 5 */
/***/ function(module, exports) {

	var Helpers = {}

	Helpers.copy = function(o) {
	  return (o instanceof Array)
	  ? o.map(Helpers.copy)
	  : (typeof o === "string" || typeof o === "number")
	  ? o
	  : copyObject(o);
	}

	function copyObject(o) {
	  var obj = {};
	  for (var k in o) obj[k] = Helpers.copy(o[k]);
	  return obj;
	}

	//Grouping cosArctan and sinArctan (see below)
	Helpers.arctans = function(dx, dy){
	  var div = dx/dy,
	      sqrt = Math.sqrt(1+(div*div)),
	      signedSqrt = (dy > 0) ? sqrt : -sqrt,
	      cos = 1 / signedSqrt,
	      sin = div * cos;

	  return {
	    cos: cos,
	    sin: sin
	  }
	}

	/*

	function cosArctan(dx,dy){
	  var div = dx/dy;
	  return (dy>0)?
	  (1/Math.sqrt(1+(div*div))):
	  (-1/Math.sqrt(1+(div*div)));
	}

	function sinArctan(dx,dy){
	  var div = dx/dy;
	  return (dy>0)?
	  (div/Math.sqrt(1+(div*div))):
	  (-div/Math.sqrt(1+(div*div)));
	}

	*/

	// TODO what's this method for ? To optimize
	Helpers.object = function(arcs, o) {
	  function arc(i, points) {
	    if (points.length) points.pop();
	    for (var a = arcs[i < 0 ? ~i : i], k = 0, n = a.length; k < n; ++k) {
	      points.push(a[k]);
	    }
	    if (i < 0) reverse(points, n);
	  }

	  function line(arcs) {
	    var points = [];
	    for (var i = 0, n = arcs.length; i < n; ++i) arc(arcs[i], points);
	    return points;
	  }

	  function polygon(arcs) {
	    return arcs.map(line);
	  }

	  function geometry(o) {
	    o = Object.create(o);
	    o.properties = o.properties; // TODO possible duplicate
	    o.coordinates = geometryType[o.type](o.arcs);
	    //type is in o's prototype, which will be lost by worker.postMessage
	    o.type = o.type
	    return o;
	  }
	  var geometryType = {
	    LineString: line,
	    MultiLineString: polygon,
	    Polygon: polygon,
	    MultiPolygon: function(arcs) { return arcs.map(polygon); }
	  };

	  return o.type === "GeometryCollection"
	  ? (o = Object.create(o), o.geometries = o.geometries.map(geometry), o)
	  : geometry(o);
	}


	function reverse(array, n) {
	  var t, j = array.length, i = j - n; while (i < --j) t = array[i], array[i++] = array[j], array[j] = t;
	}

	Helpers.properties = function(obj) {
	  return obj.properties || {};
	}

	Helpers.transformer = function(tf) {
	  var kx = tf.scale[0],
	  ky = tf.scale[1],
	  dx = tf.translate[0],
	  dy = tf.translate[1];

	  function transform(c) {
	    return [c[0] * kx + dx, c[1] * ky + dy];
	  }

	  transform.invert = function(c) {
	    return [(c[0] - dx) / kx, (c[1]- dy) / ky];
	  };

	  return transform;
	};

	module.exports = Helpers


/***/ },
/* 6 */
/***/ function(module, exports) {

	(function(root){

	  // Let's borrow a couple of things from Underscore that we'll need

	  // _.each
	  var breaker = {},
	      AP = Array.prototype,
	      OP = Object.prototype,

	      hasOwn = OP.hasOwnProperty,
	      toString = OP.toString,
	      forEach = AP.forEach,
	      indexOf = AP.indexOf,
	      slice = AP.slice;

	  var _each = function( obj, iterator, context ) {
	    var key, i, l;

	    if ( !obj ) {
	      return;
	    }
	    if ( forEach && obj.forEach === forEach ) {
	      obj.forEach( iterator, context );
	    } else if ( obj.length === +obj.length ) {
	      for ( i = 0, l = obj.length; i < l; i++ ) {
	        if ( i in obj && iterator.call( context, obj[i], i, obj ) === breaker ) {
	          return;
	        }
	      }
	    } else {
	      for ( key in obj ) {
	        if ( hasOwn.call( obj, key ) ) {
	          if ( iterator.call( context, obj[key], key, obj) === breaker ) {
	            return;
	          }
	        }
	      }
	    }
	  };

	  // _.isFunction
	  var _isFunction = function( obj ) {
	    return !!(obj && obj.constructor && obj.call && obj.apply);
	  };

	  // _.extend
	  var _extend = function( obj ) {

	    _each( slice.call( arguments, 1), function( source ) {
	      var prop;

	      for ( prop in source ) {
	        if ( source[prop] !== void 0 ) {
	          obj[ prop ] = source[ prop ];
	        }
	      }
	    });
	    return obj;
	  };

	  // $.inArray
	  var _inArray = function( elem, arr, i ) {
	    var len;

	    if ( arr ) {
	      if ( indexOf ) {
	        return indexOf.call( arr, elem, i );
	      }

	      len = arr.length;
	      i = i ? i < 0 ? Math.max( 0, len + i ) : i : 0;

	      for ( ; i < len; i++ ) {
	        // Skip accessing in sparse arrays
	        if ( i in arr && arr[ i ] === elem ) {
	          return i;
	        }
	      }
	    }

	    return -1;
	  };

	  // And some jQuery specific helpers

	  var class2type = {};

	  // Populate the class2type map
	  _each("Boolean Number String Function Array Date RegExp Object".split(" "), function(name, i) {
	    class2type[ "[object " + name + "]" ] = name.toLowerCase();
	  });

	  var _type = function( obj ) {
	    return obj == null ?
	      String( obj ) :
	      class2type[ toString.call(obj) ] || "object";
	  };

	  // Now start the jQuery-cum-Underscore implementation. Some very
	  // minor changes to the jQuery source to get this working.

	  // Internal Deferred namespace
	  var _d = {};
	  // String to Object options format cache
	  var optionsCache = {};

	  // Convert String-formatted options into Object-formatted ones and store in cache
	  function createOptions( options ) {
	    var object = optionsCache[ options ] = {};
	    _each( options.split( /\s+/ ), function( flag ) {
	      object[ flag ] = true;
	    });
	    return object;
	  }

	  _d.Callbacks = function( options ) {

	    // Convert options from String-formatted to Object-formatted if needed
	    // (we check in cache first)
	    options = typeof options === "string" ?
	      ( optionsCache[ options ] || createOptions( options ) ) :
	      _extend( {}, options );

	    var // Last fire value (for non-forgettable lists)
	      memory,
	      // Flag to know if list was already fired
	      fired,
	      // Flag to know if list is currently firing
	      firing,
	      // First callback to fire (used internally by add and fireWith)
	      firingStart,
	      // End of the loop when firing
	      firingLength,
	      // Index of currently firing callback (modified by remove if needed)
	      firingIndex,
	      // Actual callback list
	      list = [],
	      // Stack of fire calls for repeatable lists
	      stack = !options.once && [],
	      // Fire callbacks
	      fire = function( data ) {
	        memory = options.memory && data;
	        fired = true;
	        firingIndex = firingStart || 0;
	        firingStart = 0;
	        firingLength = list.length;
	        firing = true;
	        for ( ; list && firingIndex < firingLength; firingIndex++ ) {
	          if ( list[ firingIndex ].apply( data[ 0 ], data[ 1 ] ) === false && options.stopOnFalse ) {
	            memory = false; // To prevent further calls using add
	            break;
	          }
	        }
	        firing = false;
	        if ( list ) {
	          if ( stack ) {
	            if ( stack.length ) {
	              fire( stack.shift() );
	            }
	          } else if ( memory ) {
	            list = [];
	          } else {
	            self.disable();
	          }
	        }
	      },
	      // Actual Callbacks object
	      self = {
	        // Add a callback or a collection of callbacks to the list
	        add: function() {
	          if ( list ) {
	            // First, we save the current length
	            var start = list.length;
	            (function add( args ) {
	              _each( args, function( arg ) {
	                var type = _type( arg );
	                if ( type === "function" ) {
	                  if ( !options.unique || !self.has( arg ) ) {
	                    list.push( arg );
	                  }
	                } else if ( arg && arg.length && type !== "string" ) {
	                  // Inspect recursively
	                  add( arg );
	                }
	              });
	            })( arguments );
	            // Do we need to add the callbacks to the
	            // current firing batch?
	            if ( firing ) {
	              firingLength = list.length;
	            // With memory, if we're not firing then
	            // we should call right away
	            } else if ( memory ) {
	              firingStart = start;
	              fire( memory );
	            }
	          }
	          return this;
	        },
	        // Remove a callback from the list
	        remove: function() {
	          if ( list ) {
	            _each( arguments, function( arg ) {
	              var index;
	              while( ( index = _inArray( arg, list, index ) ) > -1 ) {
	                list.splice( index, 1 );
	                // Handle firing indexes
	                if ( firing ) {
	                  if ( index <= firingLength ) {
	                    firingLength--;
	                  }
	                  if ( index <= firingIndex ) {
	                    firingIndex--;
	                  }
	                }
	              }
	            });
	          }
	          return this;
	        },
	        // Control if a given callback is in the list
	        has: function( fn ) {
	          return _inArray( fn, list ) > -1;
	        },
	        // Remove all callbacks from the list
	        empty: function() {
	          list = [];
	          return this;
	        },
	        // Have the list do nothing anymore
	        disable: function() {
	          list = stack = memory = undefined;
	          return this;
	        },
	        // Is it disabled?
	        disabled: function() {
	          return !list;
	        },
	        // Lock the list in its current state
	        lock: function() {
	          stack = undefined;
	          if ( !memory ) {
	            self.disable();
	          }
	          return this;
	        },
	        // Is it locked?
	        locked: function() {
	          return !stack;
	        },
	        // Call all callbacks with the given context and arguments
	        fireWith: function( context, args ) {
	          args = args || [];
	          args = [ context, args.slice ? args.slice() : args ];
	          if ( list && ( !fired || stack ) ) {
	            if ( firing ) {
	              stack.push( args );
	            } else {
	              fire( args );
	            }
	          }
	          return this;
	        },
	        // Call all the callbacks with the given arguments
	        fire: function() {
	          self.fireWith( this, arguments );
	          return this;
	        },
	        // To know if the callbacks have already been called at least once
	        fired: function() {
	          return !!fired;
	        }
	      };

	    return self;
	  };

	  _d.Deferred = function( func ) {

	    var tuples = [
	        // action, add listener, listener list, final state
	        [ "resolve", "done", _d.Callbacks("once memory"), "resolved" ],
	        [ "reject", "fail", _d.Callbacks("once memory"), "rejected" ],
	        [ "notify", "progress", _d.Callbacks("memory") ]
	      ],
	      state = "pending",
	      promise = {
	        state: function() {
	          return state;
	        },
	        always: function() {
	          deferred.done( arguments ).fail( arguments );
	          return this;
	        },
	        then: function( /* fnDone, fnFail, fnProgress */ ) {
	          var fns = arguments;

	          return _d.Deferred(function( newDefer ) {

	            _each( tuples, function( tuple, i ) {
	              var action = tuple[ 0 ],
	                fn = fns[ i ];

	              // deferred[ done | fail | progress ] for forwarding actions to newDefer
	              deferred[ tuple[1] ]( _isFunction( fn ) ?

	                function() {
	                  var returned;
	                  try { returned = fn.apply( this, arguments ); } catch(e){
	                    newDefer.reject(e);
	                    return;
	                  }

	                  if ( returned && _isFunction( returned.promise ) ) {
	                    returned.promise()
	                      .done( newDefer.resolve )
	                      .fail( newDefer.reject )
	                      .progress( newDefer.notify );
	                  } else {
	                    newDefer[ action !== "notify" ? 'resolveWith' : action + 'With']( this === deferred ? newDefer : this, [ returned ] );
	                  }
	                } :

	                newDefer[ action ]
	              );
	            });

	            fns = null;

	          }).promise();

	        },
	        // Get a promise for this deferred
	        // If obj is provided, the promise aspect is added to the object
	        promise: function( obj ) {
	          return obj != null ? _extend( obj, promise ) : promise;
	        }
	      },
	      deferred = {};

	    // Keep pipe for back-compat
	    promise.pipe = promise.then;

	    // Add list-specific methods
	    _each( tuples, function( tuple, i ) {
	      var list = tuple[ 2 ],
	        stateString = tuple[ 3 ];

	      // promise[ done | fail | progress ] = list.add
	      promise[ tuple[1] ] = list.add;

	      // Handle state
	      if ( stateString ) {
	        list.add(function() {
	          // state = [ resolved | rejected ]
	          state = stateString;

	        // [ reject_list | resolve_list ].disable; progress_list.lock
	        }, tuples[ i ^ 1 ][ 2 ].disable, tuples[ 2 ][ 2 ].lock );
	      }

	      // deferred[ resolve | reject | notify ] = list.fire
	      deferred[ tuple[0] ] = list.fire;
	      deferred[ tuple[0] + "With" ] = list.fireWith;
	    });

	    // Make the deferred a promise
	    promise.promise( deferred );

	    // Call given func if any
	    if ( func ) {
	      func.call( deferred, deferred );
	    }

	    // All done!
	    return deferred;
	  };

	  // Deferred helper
	  _d.when = function( subordinate /* , ..., subordinateN */ ) {
	    var i = 0,
	      resolveValues = _type(subordinate) === 'array' && arguments.length === 1 ?
	        subordinate : slice.call( arguments ),
	      length = resolveValues.length,

	      // the count of uncompleted subordinates
	      remaining = length !== 1 || ( subordinate && _isFunction( subordinate.promise ) ) ? length : 0,

	      // the master Deferred. If resolveValues consist of only a single Deferred, just use that.
	      deferred = remaining === 1 ? subordinate : _d.Deferred(),

	      // Update function for both resolve and progress values
	      updateFunc = function( i, contexts, values ) {
	        return function( value ) {
	          contexts[ i ] = this;
	          values[ i ] = arguments.length > 1 ? slice.call( arguments ) : value;
	          if( values === progressValues ) {
	            deferred.notifyWith( contexts, values );
	          } else if ( !( --remaining ) ) {
	            deferred.resolveWith( contexts, values );
	          }
	        };
	      },

	      progressValues, progressContexts, resolveContexts;

	    // add listeners to Deferred subordinates; treat others as resolved
	    if ( length > 1 ) {
	      progressValues = new Array( length );
	      progressContexts = new Array( length );
	      resolveContexts = new Array( length );
	      for ( ; i < length; i++ ) {
	        if ( resolveValues[ i ] && _isFunction( resolveValues[ i ].promise ) ) {
	          resolveValues[ i ].promise()
	            .done( updateFunc( i, resolveContexts, resolveValues ) )
	            .fail( deferred.reject )
	            .progress( updateFunc( i, progressContexts, progressValues ) );
	        } else {
	          --remaining;
	        }
	      }
	    }

	    // if we're not waiting on anything, resolve the master
	    if ( !remaining ) {
	      deferred.resolveWith( resolveContexts, resolveValues );
	    }

	    return deferred.promise();
	  };

	  // Try exporting as a Common.js Module
	  if ( typeof module !== "undefined" && module.exports ) {
	    module.exports = _d;

	  // Or mixin to Underscore.js
	  } else if ( typeof root._ !== "undefined" ) {
	    root._.mixin(_d);

	  // Or assign it to window._
	  } else {
	    root._ = _d;
	  }

	})(this);


/***/ }
/******/ ])
});
;