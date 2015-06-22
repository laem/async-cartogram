What ?
-----------

A **js cartogram library** that does not block the browser, using **web workers**, forked from the original [d3 cartogram](http://prag.ma/code/d3-cartogram/). It can also be run with **node**, which is nice for preprocessing geometries offline.

- Designed to handle a series of cartograms (called tasks) : compute all the maps at once then display them smoothly.

- Using only a geo subset of d3 (you'll probably use d3 again to display the map in the browser, add colors, transitions, hovers...).

- Also attempts to add the effective range tip from Algorithms for Cartogram Animation,
by Ouyang et al., hopefully making it slightly faster.

Usage :
--------

The library is compiled in `dist/` as [UMD](https://github.com/umdjs/umd), and has been tested using `require("dist/async-cartogram.js")`

```js
var promiseOfGeometries = AsyncCartogram(
  {
    topology: topojsonData,
    // geometries to reshape:
    geometries: topojsonData.objects.OBJECTNAME.geometries,
    projection: {
      name: 'mercator',
      translation: [X,Y],
      scaling: scalingFactor, //e.g. 2000
      center: [long, lat]
    }
  },
  values, // { taskId => { geoJsonFeatureValue => newArea} }
  featureProperty // geoJsonFeatureIdKey to link geojson features to values
);
```

See the [example](https://github.com/laem/async-cartogram/tree/master/example) app (cartogram of Paris) for detailed usage.


![L'escargot déformé](https://github.com/laem/async-cartogram/tree/master/example/capture.png)

Check the project for which this lib was created, [a 90-years cartogram of europe](https://github.com/laem/eurpop).
