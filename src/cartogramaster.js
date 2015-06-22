function spawnWorker(){
  var Wok
  if (typeof window !== 'undefined'){//browser
    Wok = require('worker?inline=true!./cartogram-worker.js');
  } else {//node
    Wok = require('./cartogram-worker.js')
  }
  return new Wok()
}


var _ = require('underscore.deferred')

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
