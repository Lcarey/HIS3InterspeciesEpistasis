// some mind blowing mobile nav jquery
$(document).ready(function() {
  $("#toggle").click(function() {
    $(this).fadeToggle("slow",0);
    $(".overlay-nav").fadeIn("slow,", 0)
  $(".close").click(function() {
    $(this).fadeToggle("slow", 0)
    $(".overlay-nav").fadeOut("slow", 0)
    $("#toggle").fadeIn("slow", 0)
    $(".close").fadeIn("slow", 0)

});
});
});

// About Link going to about picture with supporting paragraphs

$(document).ready(function() {
  $('.about').click(function() {
    $('html, body').animate({
      scrollTop: $('.self-pic').offset().top
    }, 1200);
  });
});

//Portfolio Link going to flexbox portfolio section
$(document).ready(function() {
  $('.work').click(function() {
    $('html, body').animate({
      scrollTop: $('.flex').offset().top
    }, 1200);
  });
});

// Contact link going to contact section
$(document).ready(function() {
  $('.contact-nav').click(function() {
    $('html, body').animate({
      scrollTop: $('.contact').offset().top
    }, 1200);
  });
});

// Work link going to portfolio section
$(document).ready(function() {
  $('.work').click(function() {
    $('html, body').animate({
      scrollTop: $('.flex').offset().top
    }, 1200);
  });
});

$(document).ready(function() {
  $('a').click(function() {
    $('overlay-nav').fadeOut();
  });
  });

// Figure animations
$(document).ready(function() {
  $('#one').on({
    mouseenter: function(){
      $('#img1').css("height", "75px")
      $('#img1').css("opacity", "1")
    }
  });
  $('#one, #two, #three, #four, #five, #six').on({
    mouseleave: function(){
    $('figcaption').css("height", "0px")
    $('figcaption').css("opacity", "0")
    }

  });
});

$(document).ready(function() {
  $('#two').on({
    mouseenter: function(){
      $('#img2').css("height", "75px")
      $('#img2').css("opacity", "1")
    }
  })
});

$(document).ready(function() {
  $('#three').on({
    mouseenter: function(){
      $('#img3').css("height", "75px")
      $('#img3').css("opacity", "1")
    }
  })
});

$(document).ready(function() {
  $('#four').on({
    mouseenter: function(){
      $('#img4').css("height", "75px")
      $('#img4').css("opacity", "1")
    }
  })
});

$(document).ready(function() {
  $('#five').on({
    mouseenter: function(){
      $('#img5').css("height", "75px")
      $('#img5').css("opacity", "1")
    }
  })
});

$(document).ready(function() {
  $('#six').on({
    mouseenter: function(){
      $('#img6').css("height", "75px")
      $('#img6').css("opacity", "1")
    }
  })
});
// End figcaption animations

// figcaption doesn't show up under 960px. Should get better performance on mobile.
$(document).ready(function() {
  $(window).resize(function() {
 if ($(window).width() < 960) {
    $('figcaption').css('display', 'none')
 }
 else {
   $('figcaption').css('display', 'block')
}
});
});







// Give Figcaption an ID and on MouseEnter And Leave use CSS properties #work
