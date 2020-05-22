//===Update default infile values
function loadDefaultsAjax () {
  $.ajax({
  	type : "POST",
    cache: false,
  	url : "loadDefaults",
    data: {test: "test"},
  	success: function (data) {
          document.getElementById("shot").value = data["shot"];
          document.getElementById("tmin").value = data["tmin"];
          document.getElementById("tmax").value = data["tmax"];
          document.getElementById("nTrace").value = data["nTrace"];
          document.getElementById("ROIGridRes").value = data["ROIGridRes"];
          document.getElementById("gridRes").value = data["gridRes"];
          document.getElementById("STPfile").value = data["STPfile"];
          document.getElementById("lq").value = data["lq"];
          document.getElementById("S").value = data["S"];
          document.getElementById("P").value = data["Psol"];
          document.getElementById("qBG").value = data["qBG"];
          document.getElementById("xMinOF").value = data["xMin"];
          document.getElementById("xMaxOF").value = data["xMax"];
          document.getElementById("yMinOF").value = data["yMin"];
          document.getElementById("yMaxOF").value = data["yMax"];
          document.getElementById("zMinOF").value = data["zMin"];
          document.getElementById("zMaxOF").value = data["zMax"];
          document.getElementById("tMinOF").value = data["tMin"];
          document.getElementById("tMaxOF").value = data["tMax"];
          document.getElementById("deltaT").value = data["deltaT"];
          document.getElementById("writeDeltaT").value = data["writeDeltaT"];
          document.getElementById("STLscale").value = data["STLscale"];
          document.getElementById("meshMinLevel").value = data["meshMinLevel"];
          document.getElementById("meshMaxLevel").value = data["meshMaxLevel"];
          document.getElementById("xMidOF").value = data["xMid"];
          document.getElementById("yMidOF").value = data["yMid"];
          document.getElementById("zMidOF").value = data["zMid"];
          document.getElementById("xCoord_T").value = data["xProbe"];
          document.getElementById("yCoord_T").value = data["yProbe"];
          document.getElementById("zCoord_T").value = data["zProbe"];
  		},
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
  	});
};

//===Update PFC parts and intersects
function loadPFCDefaultsAjax () {
    //First set all parts' checkboxes to false
    $.ajax({
    	type : "POST",
      cache: false,
    	url : "loadPFCParts",
      data: {test: "test"},
    	success: function (data) {
        for (var i=0; i<data.length; i++){
          document.getElementById(data[i]+'1').checked = false;
          document.getElementById(data[i]+'2').checked = false;
        }
        updatePFCBoxes();
    		},
      error: function(xhr, status, error) {
        var err = eval("(" + xhr.responseText + ")");
        alert(err.Message);
        }
    	});
};
function updatePFCBoxes () {
      //Now set the correct checkboxes to true
      $.ajax({
      	type : "POST",
        cache: false,
      	url : "loadPFCDefaults",
        data: {test: "test"},
      	success: function (data) {
          var parts = data["parts"];
          var intersects = data["intersects"];

          for (var i=0; i<parts.length; i++){
            document.getElementById(parts[i]+'1').checked = true;
          }

          for (var i=0; i<intersects.length; i++){
            document.getElementById(intersects[i]+'2').checked = true;
          }

      		},
        error: function(xhr, status, error) {
          var err = eval("(" + xhr.responseText + ")");
          alert(err.Message);
          }
      	});
};

// This is the button click event
$('#loadDefaults').click(function () {
  loadDefaultsAjax();
  loadPFCDefaultsAjax();
});
