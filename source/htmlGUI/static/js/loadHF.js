function loadHFajax () {
  var e = document.getElementById("HFdropdown");
  var mode =  e.options[e.selectedIndex].value;

  if (mode == 'multiExp') {
    var lqPN = document.getElementById("lqPN2").value;
    var lqPF = document.getElementById("lqPF2").value;
    var lqCN = document.getElementById("lqCN2").value;
    var lqCF = document.getElementById("lqCF2").value;
    var fracPN = document.getElementById("fracPN2").value;
    var fracPF = document.getElementById("fracPF2").value;
    var fracCN = document.getElementById("fracCN2").value;
    var fracCF = document.getElementById("fracCF2").value;
    var p = document.getElementById("P2").value;
    var LRmask = document.getElementById("LRmask2").checked;
    var LRpower = document.getElementById("LRpower2").value;
  } else if (mode == 'limiter') {
    var lqCN = document.getElementById("lqCN3").value;
    var lqCF = document.getElementById("lqCF3").value;
    var fracCN = document.getElementById("fracCN3").value;
    var fracCF = document.getElementById("fracCF3").value;
    var p = document.getElementById("P3").value;
    var LRmask = document.getElementById("LRmask3").checked;
    var LRpower = document.getElementById("LRpower3").value;
  } else {
    var lqEich = document.getElementById("lqEich").value;
    var s = document.getElementById("S").value;
    var p = document.getElementById("P1").value;
    var qBG = document.getElementById("qBG").value;
    var LRmask = document.getElementById("LRmask1").checked;
    var LRpower = document.getElementById("LRpower1").value;
  }



  $.ajax({
    type : "POST",
    cache: false,
    url : "loadHF",
    data: {lqEich: lqEich, s: s, p: p, qBG: qBG, lqPN: lqPN, lqPF: lqPF, lqCN: lqCN, lqCF: lqCF,
                      fracPN: fracPN, fracPF: fracPF, fracCN: fracCN, fracCF: fracCF, mode: mode,
                      LRmask: LRmask, LRpower: LRpower},
    success: function (data) {
//      alert('HF Success')
      },
    error: function(xhr, status, error) {
      var err = eval("(" + xhr.responseText + ")");
      alert(err.Message);
      }
    });
}
$('#loadHF1').click(function () {
  loadHFajax();
});
$('#loadHF2').click(function () {
  loadHFajax();
});
$('#loadHF3').click(function () {
  loadHFajax();
});
