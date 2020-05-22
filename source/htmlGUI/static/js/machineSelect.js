function machineSelectAjax () {

  if (document.getElementById('nstxBox').checked) {
  machine = document.getElementById('nstxBox').value;
} else if (document.getElementById('st40Box').checked) {
  machine = document.getElementById('st40Box').value;
} else if (document.getElementById('d3dBox').checked) {
  machine = document.getElementById('d3dBox').value;
}

  $.ajax({
    type : "POST",
    cache: false,
    url : "machineSelect",
    data: {machine:machine},
    success: function (data) {
        window.location='run';      
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
    });
}
$('#machineSelect').click(function () {
  machineSelectAjax();
});
