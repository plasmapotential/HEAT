function loadPFCParts() {
  $.ajax({
    type : "POST",
    cache: false,
    url : "loadPFCParts",
    success: function (data) {
      var div1 = document.getElementById("scrolldiv1");
      var div2 = document.getElementById("scrolldiv2");
      for (var i=0; i<data.length; i++){
        //These next two lines are ugly but they work
        //We dynamically create a checkbox list for HF parts and intersect parts
        //that we get out of the partlist file.  We call a flask function that
        //returns a list of these parts (as strings), then dynamically insert them
        //into the div1 and div2 <div>s as raw html.
        div1.innerHTML += '<div class="field"><input type="checkbox" id="'+data[i]+'1"><label for="'+data[i]+'1">'+data[i]+'</label></div>';
        div2.innerHTML += '<div class="field"><input type="checkbox" id="'+data[i]+'2"><label for="'+data[i]+'2">'+data[i]+'</label></div>';

      }
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
    }
    });
};

// This is the button click event
$( document ).ready(function () {
  loadPFCParts();
});
