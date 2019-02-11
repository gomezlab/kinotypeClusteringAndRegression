(window["webpackJsonp"] = window["webpackJsonp"] || []).push([[1],{

/***/ "./node_modules/sigma/build/plugins/sigma.parsers.json.min.js":
/*!********************************************************************!*\
  !*** ./node_modules/sigma/build/plugins/sigma.parsers.json.min.js ***!
  \********************************************************************/
/*! no static exports found */
/***/ (function(module, exports) {

eval("/*** IMPORTS FROM imports-loader ***/\n(function() {\n\n(function(a){\"use strict\";if(\"undefined\"==typeof sigma)throw\"sigma is not declared\";sigma.utils.pkg(\"sigma.parsers\"),sigma.utils.pkg(\"sigma.utils\"),sigma.utils.xhr=function(){if(window.XMLHttpRequest)return new XMLHttpRequest;var a,b;if(window.ActiveXObject){a=[\"Msxml2.XMLHTTP.6.0\",\"Msxml2.XMLHTTP.3.0\",\"Msxml2.XMLHTTP\",\"Microsoft.XMLHTTP\"];for(b in a)try{return new ActiveXObject(a[b])}catch(a){}}return null},sigma.parsers.json=function(a,b,c){var d,e=sigma.utils.xhr();if(!e)throw\"XMLHttpRequest not supported, cannot load the file.\";e.open(\"GET\",a,!0),e.onreadystatechange=function(){4===e.readyState&&(d=JSON.parse(e.responseText),b instanceof sigma?(b.graph.clear(),b.graph.read(d)):\"object\"==typeof b?(b.graph=d,b=new sigma(b)):\"function\"==typeof b&&(c=b,b=null),c&&c(b||d))},e.send()}}).call(this);\n}.call(window));\n\n//# sourceURL=webpack:///./node_modules/sigma/build/plugins/sigma.parsers.json.min.js?");

/***/ })

}]);