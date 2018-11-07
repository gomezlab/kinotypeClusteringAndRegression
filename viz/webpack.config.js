const path = require('path');

module.exports = {
  mode: 'development',
  entry: './src/index.js',
  output: {
    filename: 'main.js',
    path: path.resolve(__dirname, 'dist')
  },
  resolve: {
    modules: [
      'node_modules',
    ],
  },
  module: {
    rules: [
      {
        test: /sigma.*/,
        use: 'imports-loader?this=>window',
      },
    ],
  }
};
