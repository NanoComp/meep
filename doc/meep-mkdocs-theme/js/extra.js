//HR 20170807 startup javascript code to fix broken search box on readthedocs
//inspired by similar code in the nodecmu project:
//https://github.com/nodemcu/nodemcu-firmware/blob/master/docs/js/extra.js

  /**
   * Analyzes the URL of the current page to find out what the selected GitHub branch is. It's usually
   * part of the location path. The code needs to distinguish between running MkDocs standalone
   * and docs served from RTD. If no valid branch could be determined 'dev' returned.
   *
   * @returns GitHub branch name
   */
  function determineSelectedBranch() {
    var branch = 'dev', path = window.location.pathname;
    if (window.location.origin.indexOf('readthedocs') > -1) {
      // path is like /en/<branch>/<lang>/build/ -> extract 'lang'
      // split[0] is an '' because the path starts with the separator
      var thirdPathSegment = path.split('/')[2];
      // 'latest' is an alias on RTD for the 'dev' branch - which is the default for 'branch' here
      if (thirdPathSegment != 'latest') {
        branch = thirdPathSegment;
      }
    }
    return branch;
  }

 function fixSearch() {
    var target = document.getElementById('rtd-search-form');
    var config = {attributes: true, childList: true};

    var observer = new MutationObserver(function(mutations) {
      // if it isn't disconnected it'll loop infinitely because the observed element is modified
      observer.disconnect();
      var form = $('#rtd-search-form');
      form.empty();

 // HR 20170807 the following fancier code attempts to determine the appropriate github branch/readthedocs version on the fly; for now I'm hard-coding it to 'latest'
 //     form.attr('action', 'https://' + window.location.hostname + '/en/' + determineSelectedBranch() + '/search.html');
        form.attr('action', 'https://' + window.location.hostname + '/en/latest/search.html');

      $('<input>').attr({
        type: "text",
        name: "q",
        placeholder: "Search docs"
      }).appendTo(form);
    });
    // don't run this outside RTD hosting
    if (window.location.origin.indexOf('readthedocs') > -1) {
      observer.observe(target, config);
    }
  }

$( document ).ready(function() {
    fixSearch();
});
