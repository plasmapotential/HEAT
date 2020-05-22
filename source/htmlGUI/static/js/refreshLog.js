function refreshLogAjax() {
  $.ajax({
    type : "POST",
    cache: false,
    url : "refreshLog",
    success: function (logData) {
      document.getElementById("log").innerHTML = logData;
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
    }
    });
};

// This is the button click event
$( document ).ready(function () {
  setInterval(refreshLogAjax, 3000)
});
