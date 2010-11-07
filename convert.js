var fs = require('fs'),
    sys = require('sys');
var f = fs.readFileSync('potrace.js', 'utf-8');
sys.print(f.split(';').map(function(line){
  //sys.print(line);
  
  if(/\w+ \w+\(.*?\)\s+\{/.test(line)){
    //console.log('function declaration');
    //console.log(line);
    line = line.replace(/(\w+) (\w+)\((.*?)\)(\s+\{)/g, function(all, type, name, args, trail){
      return 'function /*'+type+'*/ '+name+'('+args.split(',').map(function(a){
        return a.replace(/(\w+) (\w+)/g, '/*$1*/ $2')
      }).join(',')+')'+trail
    })
    //console.log(line);
  }
  if(/\w+ \w+ =/.test(line)){
    //console.log('declaration');
    //console.log(line);
    line = line.replace(/(\w+) (\w+ =)/g,'var $2')
  }
  if(/\w+ \w+, \w+/.test(line)){
    //console.log('declaration');
    //console.log(line);
    line = line.replace(/(\w+) (\w+, \w+)/g,'var $2')
  }
  return line
}).join(';').split('\n').filter(function(a){
  return !/\/\/\//.test(a);
}).join('\n'))
