function loadCADajax() {
  var ROIGridRes = document.getElementById("ROIGridRes").value;
  var gridRes = document.getElementById("gridRes").value;
  var STPfile = document.getElementById("STPfile").value;
  $.ajax({
    type : "POST",
    cache: false,
    url : "loadCAD",
    data: {ROIGridRes: ROIGridRes, gridRes: gridRes, STPfile: STPfile},
    success: function (data) {
//      alert('CAD SUCCESS')
//      $("#eq").attr('src',"data:image/png;base64," + png);
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
    }
    });
};

// This is the button click event
$('#loadCAD').click(function () {
  loadCADajax();
});
