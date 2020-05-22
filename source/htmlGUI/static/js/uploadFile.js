function loadFile () {
  var file = document.getElementById('fileChoose').files[0];
  var reader = new FileReader();
  reader.readAsText(file, 'UTF-8');
  reader.onload = shipOff;
}

function shipOff(event){
  var result = event.target.result;
  var fileName = document.getElementById('fileChoose').files[0].name;
  $.ajax({
    type : "POST",
    cache: false,
    url : "uploadFile",
    data: {result: result, name: fileName},
    success: function (data) {
      document.getElementById("shot").value = data["shot"];
      document.getElementById("tmin").value = data["tmin"];
      document.getElementById("tmax").value = data["tmax"];
      document.getElementById("ROIGridRes").value = data["ROIGridRes"];
      document.getElementById("gridRes").value = data["gridRes"];
      document.getElementById("STPfile").value = data["STPfile"];
      document.getElementById("lq").value = data["lq"];
      document.getElementById("S").value = data["S"];
      document.getElementById("P").value = data["Psol"];
      document.getElementById("qBG").value = data["qBG"];
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
    });
};

// This is the button click event
$('#uploadFile').click(function () {
  loadFile();
});
