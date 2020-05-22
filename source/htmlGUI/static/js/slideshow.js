function openSlide(evt, slideName) {
  // Declare all variables
  var i, slides, dots;

  // Get all elements with class="tabcontent" and hide them
  slides = document.getElementsByClassName("slide-container");
  for (i = 0; i < slides.length; i++) {
    slides[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  dots = document.getElementsByClassName("dot");
  for (i = 0; i < dots.length; i++) {
    dots[i].className = dots[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(slideName).style.display = "block";
  evt.currentTarget.className += " active";
}
