function loadMHDajax () {
  var shot = document.getElementById("shot").value;
  var tmin = document.getElementById("tmin").value;
  var tmax = document.getElementById("tmax").value;
  var nTrace = document.getElementById("nTrace").value;
  var ionDirection = document.getElementById("ionDirection").value;
  var gfilePath = document.getElementById("gfilePath").value;
  var plasma3Dmask = document.getElementById("plasma3Dmask").checked;
//  var docWidth = $(document).width();
//  var width = docWidth * 0.65 * 0.55; //gridcontainer1 * gridcontainer2
//  var height = document.getElementById("item10").offsetHeight;
//  var element = document.getElementById('item10');
//  var positionInfo = element.getBoundingClientRect();
//  var height = positionInfo.height;
//  var width = positionInfo.width;

  var height = $("#itemMHD1").height();

  $.ajax({
    type : "POST",
    cache: false,
    url : "loadMHD",
    data: {shot: shot, tmin: tmin, tmax: tmax, height: height, nTrace: nTrace,
           ionDirection: ionDirection, gfilePath: gfilePath, plasma3Dmask: plasma3Dmask},
    success: function (data) {
    //  $("#eq").attr('src',"data:image/png;base64," + png);
    //$("#eq").attr('data', "{{ url_for('static',filename='html/EQ2D.html') }}");
//    var obj = document.getElementById("eq22");
//    obj.setAttribute('data', '/u/tlooby/source/HEAT/rev6/htmlGUI/templates/EQ2D.html');
//    var cl = obj.cloneNode(true);
//    var parent = obj.parentNode;
//    parent.removeChild(obj);
//    parent.appendChild(cl);
    document.getElementById('eq').contentDocument.location.reload(true);
    document.getElementById('eq2').contentDocument.location.reload(true);
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
    });
};

// This is the button click event
$('#loadMHD').click(function () {
 loadMHDajax();
});
