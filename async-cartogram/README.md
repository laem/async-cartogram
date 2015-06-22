A js cartogram that does not block the browser using web workers, forked from the original http://prag.ma/code/d3-cartogram/.
You can also run it with nodejs.

- Designed to compute a series of cartograms (called tasks), ideally to present maps after distortions have been computed.

- For this, only a subset of d3 is used, and the usage is more constrained.

- Also attempts to add the effective range tip from Algorithms for Cartogram Animation,
by Ouyang et al., which should make it slightly faster than the original.

Usage :
--------

var promiseOfGeos = cartogramaster(
  {
    topology: topojsonData,
    // the geometries of the GeometryCollection object to reshape :
    geometries: topojsonData.objects.OBJECTNAME.geometries,
    projection: {
      name: 'mercator',
      translation: [X,Y],
      scaling: scaling
    }
  },
  values, // { taskId => { geoJsonFeatureValue => area wanted} }
  featureProperty // geoJsonFeatureIdKey to link geojson features to values
);
