function loadHFajax () {
  var lq = document.getElementById("lq").value;
  var s = document.getElementById("S").value;
  var p = document.getElementById("P").value;
  var qBG = document.getElementById("qBG").value;
  $.ajax({
    type : "POST",
    cache: false,
    url : "loadHF",
    data: {lq: lq, s: s, p: p, qBG: qBG},
    success: function (data) {
//      alert('HF Success')
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
    });
}
$('#loadHF').click(function () {
  loadHFajax();
});
