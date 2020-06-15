$('#HFdropdown').change(function () {
  var e = document.getElementById("HFdropdown");
  var prof =  e.options[e.selectedIndex].value;
  if(prof == "limiter") {
      document.getElementById('eich').style.display = 'none';
      document.getElementById('multiExp').style.display = 'none';
      document.getElementById('limiter').style.display = 'block';
  } else if(prof == 'multiExp') {
    document.getElementById('eich').style.display = 'none';
    document.getElementById('multiExp').style.display = 'block';
    document.getElementById('limiter').style.display = 'none';
  } else if(prof == 'eich') {
    document.getElementById('eich').style.display = 'block';
    document.getElementById('multiExp').style.display = 'none';
    document.getElementById('limiter').style.display = 'none';
  } else {
    document.getElementById('eich').style.display = 'none';
    document.getElementById('multiExp').style.display = 'none';
    document.getElementById('limiter').style.display = 'none';
  }
});

//onload
$(function(){
  document.getElementById('eich').style.display = 'none';
  document.getElementById('multiExp').style.display = 'none';
  document.getElementById('limiter').style.display = 'none';

});
