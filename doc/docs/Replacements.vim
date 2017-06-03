:%s/[(]\/w:\([^ ]*\) "wikilink"/(https:\/\/en.wikipedia.org\/wiki\/\1/g
:%s/\[\(.*\)\][(]Image:\([^ ]*\)[)]/\r![\1](images\/\2)\r\r/g
:%s/\[\(.*\)\][(]\/Image:\([^ ]*\) "wikilink"[)]/\r![\1](images\/\2)\r\r/g
:%s/:Category://g
:%s/(\/[Ll]ibctl "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/Libctl)/g
:%s/(\/[Ll]ibctl_manual "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/Libctl_manual)/g
:%s/(\/[Ll]ibctl_User_Reference "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/Libctl_User_Reference)/g
:%s/(\/[Ll]ibctl_tutorial "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/Libctl_tutorial)/g
:%s/(\/MPB "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/MPB)/g
:%s/(\/MPB_Installation "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/MPB_Installation)/g
:%s/(\/MPB_manual "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/MPB_manual)/g
:%s/(\/[Hh]arminv "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/harminv)/g
:%s/(\/[Hh]arminv_installation "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/Harminv_installation)/g
:%s/(\/[Hh]5utils "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/H5utils)/g
:%s/(\/Guile_and_Scheme_links "wikilink")/(http:\/\/ab-initio.mit.edu\/wiki\/index.php\/Guile_and_Scheme_links)/g
:%s/perfectly_matched_layers/Perfectly_matched_layers/g
:%s/perfectly_matched_layer/Perfectly_matched_layer/g
:%s/[(][/]Category:\([^ ]*\) "wikilink"/(\1.md/g
:%s/[(][/]\([^ ]*\) "wikilink"/(\1.md/g
:%s/\\\[/$$/g
:%s/\\\]/$$/g
:%s/\\(/$/g
:%s/\\)/$/g
:%s/^`\(.*\)`$/```\r\1\r```\r/g
:%s/\\mathbf/\\textbf/g
:%s/```\n\n```\n//g
