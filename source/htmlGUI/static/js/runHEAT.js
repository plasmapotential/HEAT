function runHEATajax() {
  var Bpc = document.getElementById("Bpc").checked;
  var Normpc = document.getElementById("Normpc").checked;
  var shadow = document.getElementById("shadow").checked;
  var HFpc = document.getElementById("HFpc").checked;
  var Btrace = document.getElementById("Btrace").checked;
  var psiPC = document.getElementById("psiPC").checked;
  var XCoord = document.getElementById("xCoord").value;
  var YCoord = document.getElementById("yCoord").value;
  var ZCoord = document.getElementById("zCoord").value;
  //openfoam variables
  var Thermal = document.getElementById("Thermal").checked;
  var Tprobe = document.getElementById("Tprobe").checked;
  var xProbe = document.getElementById("xProbe").value;
  var yProbe = document.getElementById("yProbe").value;
  var zProbe = document.getElementById("zProbe").value;
  $.ajax({
  	type : "POST",
    cache: false,
  	url : "runHEAT",
    data: {Bpc: Bpc, Btrace: Btrace, Normpc: Normpc, shadow: shadow,
          HFpc:HFpc, psiPC: psiPC,
          XCoord:XCoord, YCoord:YCoord, ZCoord:ZCoord,
          thermal: Thermal, Tprobe: Tprobe, xProbe:xProbe, yProbe:yProbe, zProbe:zProbe},
  	success: function (data) {
      alert("HEAT Run Complete");
  		},
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
  	});
};

// This is the button click event
$('#runHEAT').click(function () {
  runHEATajax();
});
