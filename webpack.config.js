var webpack = require('webpack');

module.exports = {
  entry: "./src/cartogramaster.js",
  output: {
    path: "./dist",
    filename: "async-cartogram.js",
    libraryTarget: 'umd',
    library: "AsyncCartogram"
  },

  plugins: [
    new webpack.optimize.UglifyJsPlugin({
      compress: {
        warnings: false
      }
    })
  ]

};
