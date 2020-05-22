function getAllPFCsandWrite() {
  //First get list of all the parts that we have dynamically included in checkbox areas
  $.ajax({
    type : "POST",
    cache: false,
    url : "loadPFCParts",
    data: {test: "test"},
    success: function (data) {
      var parts = [];
      var intersects = [];
      for (var i=0; i<data.length; i++){
        if (document.getElementById(data[i]+'1').checked){
          parts += data[i];
        }
        if (document.getElementById(data[i]+'2').checked){
          intersects += data[i];
        }
      }
      var dict = {parts: parts, intersects: intersects};
      writePFCsAjax(dict);
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
    });
};
function writePFCsAjax(dict) {
  //Now write PartsFile and IntersectFile based upon whats checked in GUI
  $.ajax({
    type : "POST",
    cache: false,
    url : "writePFCs",
    data: dict,
    success: function (data) {
      //alert('Wrote PFCs');
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
    }
    });
};

// This is the button click event
$('#writePFCs').click(function () {
  getAllPFCsandWrite();
});
