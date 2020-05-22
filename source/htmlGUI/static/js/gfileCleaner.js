function loadGfileAjax () {
  var gfilePath = document.getElementById("gfilePath1").value;
  var height = $("#item12").height();
  $.ajax({
    type : "POST",
    cache: false,
    url : "loadMHD",
    data: {gfilePath: gfilePath, height: height},
    success: function (png) {
      document.getElementById('eq').contentDocument.location.reload(true);
      document.getElementById('eq2').contentDocument.location.reload(true);
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
    });
}

function gfileCleanAjax () {
  var height = $("#item12").height();
  var psiRZ = document.getElementById("psiRZ").value;
  var psiSep = document.getElementById("psiSep").value;
  var psiAxis = document.getElementById("psiAxis").value;
  var Fpol = document.getElementById("Fpol").value;

  $.ajax({
    type : "POST",
    cache: false,
    url : "gfileClean",
    data: {height: height, psiRZ: psiRZ, psiSep: psiSep, psiAxis: psiAxis, Fpol: Fpol},
    success: function (png) {
      $("#eqClean").attr('src',"data:image/png;base64," + png);
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
    });
}

function writeGfileAjax () {
  var newShot = document.getElementById("newShot").value;
  var newTime = document.getElementById("newTime").value;
  var gfilePath1New = document.getElementById("gfilePath1New").value;

  $.ajax({
    type : "POST",
    cache: false,
    url : "writeGfile",
    data: {newShot: newShot, newTime: newTime, newGfile: gfilePath1New},
    success: function (data) {
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
    });
}






$('#loadGfile').click(function () {
  loadGfileAjax();
});

$('#gfileClean').click(function () {
  gfileCleanAjax();
});

$('#writeGfile').click(function () {
  writeGfileAjax();
});
