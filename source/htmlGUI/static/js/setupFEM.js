function setupFEMajax() {
  var xMinOF = document.getElementById("xMinOF").value;
  var xMaxOF = document.getElementById("xMaxOF").value;
  var yMinOF = document.getElementById("yMinOF").value;
  var yMaxOF = document.getElementById("yMaxOF").value;
  var zMinOF = document.getElementById("zMinOF").value;
  var zMaxOF = document.getElementById("zMaxOF").value;
  var tMinOF = document.getElementById("tMinOF").value;
  var tMaxOF = document.getElementById("tMaxOF").value;
  var deltaT = document.getElementById("deltaT").value;
  var writeDeltaT = document.getElementById("writeDeltaT").value;
  var STLscale = document.getElementById("STLscale").value;
  var meshMinLevel = document.getElementById("meshMinLevel").value;
  var meshMaxLevel = document.getElementById("meshMaxLevel").value;
  var xMidOF = document.getElementById("xMidOF").value;
  var yMidOF = document.getElementById("yMidOF").value;
  var zMidOF = document.getElementById("zMidOF").value;
  $.ajax({
  	type : "POST",
    cache: false,
  	url : "setupOpenFOAM",
    data: {xMinOF:xMinOF,xMaxOF:xMaxOF,yMinOF:yMinOF,yMaxOF:yMaxOF,
           zMinOF:zMinOF,zMaxOF:zMaxOF,tMinOF:tMinOF,tMaxOF:tMaxOF,
           deltaT:deltaT,writeDeltaT:writeDeltaT,STLscale:STLscale,
           meshMinLevel:meshMinLevel,meshMaxLevel:meshMaxLevel,
           xMidOF:xMidOF,yMidOF:yMidOF,zMidOF:zMidOF},
  	success: function (data) {
      //alert("OpenFOAM Setup Complete");
  		},
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
  	});
};

// This is the button click event
$('#setupFEM').click(function () {
  setupFEMajax();
});
