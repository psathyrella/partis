


<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# githubog: http://ogp.me/ns/fb/githubog#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <title>klib/ksort.h at master · attractivechaos/klib · GitHub</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <link rel="logo" type="image/svg" href="https://github-media-downloads.s3.amazonaws.com/github-logo.svg" />
    <meta property="og:image" content="https://github.global.ssl.fastly.net/images/modules/logos_page/Octocat.png">
    <meta name="hostname" content="fe2.rs.github.com">
    <meta name="ruby" content="ruby 2.0.0p247-github4 (2013-06-27) [x86_64-linux]">
    <link rel="assets" href="https://github.global.ssl.fastly.net/">
    <link rel="xhr-socket" href="/_sockets" />
    
    


    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="github" name="octolytics-app-id" />

    
    
    <link rel="icon" type="image/x-icon" href="/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="pbBE8rN4Sx1x2VNA3UK3dERocvlHQhsPGbu82N1sygI=" name="csrf-token" />

    <link href="https://github.global.ssl.fastly.net/assets/github-f226abc7983f8566b17d24236adae64ba647ffea.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://github.global.ssl.fastly.net/assets/github2-4579bfc1145b0b607e0f375917bf8fb5d21cfca2.css" media="all" rel="stylesheet" type="text/css" />
    


      <script src="https://github.global.ssl.fastly.net/assets/frameworks-308994eebca96c869a2e26e031d5eecb2f9ea76c.js" type="text/javascript"></script>
      <script src="https://github.global.ssl.fastly.net/assets/github-2211ea2c22a430934fcd6a6f5360e5324b92155d.js" type="text/javascript"></script>
      
      <meta http-equiv="x-pjax-version" content="09d83e73c5793c438edbb317733aab38">

        <link data-pjax-transient rel='permalink' href='/attractivechaos/klib/blob/57dc7b07c9315ead4c106405c2e0968d8bbe2950/ksort.h'>
  <meta property="og:title" content="klib"/>
  <meta property="og:type" content="githubog:gitrepository"/>
  <meta property="og:url" content="https://github.com/attractivechaos/klib"/>
  <meta property="og:image" content="https://github.global.ssl.fastly.net/images/gravatars/gravatar-user-420.png"/>
  <meta property="og:site_name" content="GitHub"/>
  <meta property="og:description" content="klib - A standalone and lightweight C library"/>

  <meta name="description" content="klib - A standalone and lightweight C library" />

  <meta content="563093" name="octolytics-dimension-user_id" /><meta content="attractivechaos" name="octolytics-dimension-user_login" /><meta content="1251393" name="octolytics-dimension-repository_id" /><meta content="attractivechaos/klib" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="1251393" name="octolytics-dimension-repository_network_root_id" /><meta content="attractivechaos/klib" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/attractivechaos/klib/commits/master.atom" rel="alternate" title="Recent Commits to klib:master" type="application/atom+xml" />

  </head>


  <body class="logged_out page-blob  vis-public env-production ">

    <div class="wrapper">
      
      
      


      
      <div class="header header-logged-out">
  <div class="container clearfix">

    <a class="header-logo-wordmark" href="https://github.com/">
      <span class="mega-octicon octicon-logo-github"></span>
    </a>

    <div class="header-actions">
        <a class="button primary" href="/signup">Sign up</a>
      <a class="button" href="/login?return_to=%2Fattractivechaos%2Fklib%2Fblob%2Fmaster%2Fksort.h">Sign in</a>
    </div>

    <div class="command-bar js-command-bar  in-repository">


      <ul class="top-nav">
          <li class="explore"><a href="/explore">Explore</a></li>
        <li class="features"><a href="/features">Features</a></li>
          <li class="enterprise"><a href="https://enterprise.github.com/">Enterprise</a></li>
          <li class="blog"><a href="/blog">Blog</a></li>
      </ul>
        <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<input type="text" data-hotkey="/ s" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    
      data-repo="attractivechaos/klib"
      data-branch="master"
      data-sha="8f3432a9764171b904904821c570b90d56f10469"
  >

    <input type="hidden" name="nwo" value="attractivechaos/klib" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item js-this-repository-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item js-all-repositories-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="octicon help tooltipped downwards" title="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
    </div>

  </div>
</div>


      


          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        

<ul class="pagehead-actions">


  <li>
  <a href="/login?return_to=%2Fattractivechaos%2Fklib"
  class="minibutton with-count js-toggler-target star-button entice tooltipped upwards"
  title="You must be signed in to use this feature" rel="nofollow">
  <span class="octicon octicon-star"></span>Star
</a>
<a class="social-count js-social-count" href="/attractivechaos/klib/stargazers">
  376
</a>

  </li>

    <li>
      <a href="/login?return_to=%2Fattractivechaos%2Fklib"
        class="minibutton with-count js-toggler-target fork-button entice tooltipped upwards"
        title="You must be signed in to fork a repository" rel="nofollow">
        <span class="octicon octicon-git-branch"></span>Fork
      </a>
      <a href="/attractivechaos/klib/network" class="social-count">
        43
      </a>
    </li>
</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo"></span>
          <span class="author">
            <a href="/attractivechaos" class="url fn" itemprop="url" rel="author"><span itemprop="title">attractivechaos</span></a></span
          ><span class="repohead-name-divider">/</span><strong
          ><a href="/attractivechaos/klib" class="js-current-repository js-repo-home-link">klib</a></strong>

          <span class="page-context-loader">
            <img alt="Octocat-spinner-32" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">

      <div class="repository-with-sidebar repo-container ">

        <div class="repository-sidebar">
            

<div class="repo-nav repo-nav-full js-repository-container-pjax js-octicon-loaders">
  <div class="repo-nav-contents">
    <ul class="repo-menu">
      <li class="tooltipped leftwards" title="Code">
        <a href="/attractivechaos/klib" aria-label="Code" class="js-selected-navigation-item selected" data-gotokey="c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_tags repo_branches /attractivechaos/klib">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

        <li class="tooltipped leftwards" title="Issues">
          <a href="/attractivechaos/klib/issues" aria-label="Issues" class="js-selected-navigation-item js-disable-pjax" data-gotokey="i" data-selected-links="repo_issues /attractivechaos/klib/issues">
            <span class="octicon octicon-issue-opened"></span> <span class="full-word">Issues</span>
            <span class='counter'>5</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>

      <li class="tooltipped leftwards" title="Pull Requests"><a href="/attractivechaos/klib/pulls" aria-label="Pull Requests" class="js-selected-navigation-item js-disable-pjax" data-gotokey="p" data-selected-links="repo_pulls /attractivechaos/klib/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>0</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped leftwards" title="Wiki">
          <a href="/attractivechaos/klib/wiki" aria-label="Wiki" class="js-selected-navigation-item " data-pjax="true" data-selected-links="repo_wiki /attractivechaos/klib/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
    </ul>
    <div class="repo-menu-separator"></div>
    <ul class="repo-menu">

      <li class="tooltipped leftwards" title="Pulse">
        <a href="/attractivechaos/klib/pulse" aria-label="Pulse" class="js-selected-navigation-item " data-pjax="true" data-selected-links="pulse /attractivechaos/klib/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Graphs">
        <a href="/attractivechaos/klib/graphs" aria-label="Graphs" class="js-selected-navigation-item " data-pjax="true" data-selected-links="repo_graphs repo_contributors /attractivechaos/klib/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Network">
        <a href="/attractivechaos/klib/network" aria-label="Network" class="js-selected-navigation-item js-disable-pjax" data-selected-links="repo_network /attractivechaos/klib/network">
          <span class="octicon octicon-git-branch"></span> <span class="full-word">Network</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

    </ul>

  </div>
</div>

            <div class="only-with-full-nav">
              

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=clone">
  <h3><strong>HTTPS</strong> clone URL</h3>

  <input type="text" class="clone js-url-field"
         value="https://github.com/attractivechaos/klib.git" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/attractivechaos/klib.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=clone">
  <h3><strong>Subversion</strong> checkout URL</h3>

  <input type="text" class="clone js-url-field"
         value="https://github.com/attractivechaos/klib" readonly="readonly">

  <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/attractivechaos/klib" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
</div>



<p class="clone-options">You can clone with
    <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
    <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>,
  and <a href="https://help.github.com/articles/which-remote-url-should-i-use">other methods.</a>
</p>



                <a href="/attractivechaos/klib/archive/master.zip"
                   class="minibutton sidebar-button"
                   title="Download this repository as a zip file"
                   rel="nofollow">
                  <span class="octicon octicon-cloud-download"></span>
                  Download ZIP
                </a>
            </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          


<!-- blob contrib key: blob_contributors:v21:9e61a653fc3c7e0c0e3868e5bf48e5cd -->
<!-- blob contrib frag key: views10/v8/blob_contributors:v21:9e61a653fc3c7e0c0e3868e5bf48e5cd -->

<p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

<a href="/attractivechaos/klib/find/master" data-pjax data-hotkey="t" style="display:none">Show File Finder</a>

<div class="file-navigation">
  


<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target" data-hotkey="w"
    data-master-branch="master"
    data-ref="master">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax>

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-remove-close js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/attractivechaos/klib/blob/ksw-reduce8/ksort.h" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="ksw-reduce8" data-skip-pjax="true" rel="nofollow" title="ksw-reduce8">ksw-reduce8</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/attractivechaos/klib/blob/master/ksort.h" class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target" data-name="master" data-skip-pjax="true" rel="nofollow" title="master">master</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/attractivechaos/klib" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">klib</span></a></span></span><span class="separator"> / </span><strong class="final-path">ksort.h</strong> <span class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="ksort.h" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>


  
  <div class="commit file-history-tease">
    <img class="main-avatar" height="24" src="https://1.gravatar.com/avatar/9bad6dc5cd14c6e9326a81d91b31eabb?d=https%3A%2F%2Fa248.e.akamai.net%2Fassets.github.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png&amp;s=140" width="24" />
    <span class="author"><span rel="author">Heng Li</span></span>
    <time class="js-relative-date" datetime="2011-04-10T20:29:45-07:00" title="2011-04-10 20:29:45">April 10, 2011</time>
    <div class="commit-title">
        <a href="/attractivechaos/klib/commit/5b2ce9cd462b9881c304844bc879cbef847d9e28" class="message" data-pjax="true" title="Added ks_shuffle() and ks_sample()">Added ks_shuffle() and ks_sample()</a>
    </div>

    <div class="participation">
      <p class="quickstat"><a href="#blob_contributors_box" rel="facebox"><strong>1</strong> contributor</a></p>
      
    </div>
    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list">
        <li class="facebox-user-list-item">
          <img height="24" src="https://2.gravatar.com/avatar/a750242db0e06a602a55054987a74f6d?d=https%3A%2F%2Fidenticons.github.com%2F586959b9f76621737e05744b6e5a32de.png&amp;s=140" width="24" />
          <a href="/attractivechaos">attractivechaos</a>
        </li>
      </ul>
    </div>
  </div>


<div id="files" class="bubble">
  <div class="file">
    <div class="meta">
      <div class="info">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">file</span>
          <span>299 lines (267 sloc)</span>
        <span>10.517 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
              <a class="minibutton disabled js-entice" href=""
                 data-entice="You must be signed in to make or propose changes">Edit</a>
          <a href="/attractivechaos/klib/raw/master/ksort.h" class="button minibutton " id="raw-url">Raw</a>
            <a href="/attractivechaos/klib/blame/master/ksort.h" class="button minibutton ">Blame</a>
          <a href="/attractivechaos/klib/commits/master/ksort.h" class="button minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->
            <a class="minibutton danger empty-icon js-entice" href=""
               data-entice="You must be signed in and on a branch to make or propose changes">
            Delete
          </a>
      </div><!-- /.actions -->

    </div>
        <div class="blob-wrapper data type-c js-blob-data">
        <table class="file-code file-diff">
          <tr class="file-code-line">
            <td class="blob-line-nums">
              <span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>
<span id="L128" rel="#L128">128</span>
<span id="L129" rel="#L129">129</span>
<span id="L130" rel="#L130">130</span>
<span id="L131" rel="#L131">131</span>
<span id="L132" rel="#L132">132</span>
<span id="L133" rel="#L133">133</span>
<span id="L134" rel="#L134">134</span>
<span id="L135" rel="#L135">135</span>
<span id="L136" rel="#L136">136</span>
<span id="L137" rel="#L137">137</span>
<span id="L138" rel="#L138">138</span>
<span id="L139" rel="#L139">139</span>
<span id="L140" rel="#L140">140</span>
<span id="L141" rel="#L141">141</span>
<span id="L142" rel="#L142">142</span>
<span id="L143" rel="#L143">143</span>
<span id="L144" rel="#L144">144</span>
<span id="L145" rel="#L145">145</span>
<span id="L146" rel="#L146">146</span>
<span id="L147" rel="#L147">147</span>
<span id="L148" rel="#L148">148</span>
<span id="L149" rel="#L149">149</span>
<span id="L150" rel="#L150">150</span>
<span id="L151" rel="#L151">151</span>
<span id="L152" rel="#L152">152</span>
<span id="L153" rel="#L153">153</span>
<span id="L154" rel="#L154">154</span>
<span id="L155" rel="#L155">155</span>
<span id="L156" rel="#L156">156</span>
<span id="L157" rel="#L157">157</span>
<span id="L158" rel="#L158">158</span>
<span id="L159" rel="#L159">159</span>
<span id="L160" rel="#L160">160</span>
<span id="L161" rel="#L161">161</span>
<span id="L162" rel="#L162">162</span>
<span id="L163" rel="#L163">163</span>
<span id="L164" rel="#L164">164</span>
<span id="L165" rel="#L165">165</span>
<span id="L166" rel="#L166">166</span>
<span id="L167" rel="#L167">167</span>
<span id="L168" rel="#L168">168</span>
<span id="L169" rel="#L169">169</span>
<span id="L170" rel="#L170">170</span>
<span id="L171" rel="#L171">171</span>
<span id="L172" rel="#L172">172</span>
<span id="L173" rel="#L173">173</span>
<span id="L174" rel="#L174">174</span>
<span id="L175" rel="#L175">175</span>
<span id="L176" rel="#L176">176</span>
<span id="L177" rel="#L177">177</span>
<span id="L178" rel="#L178">178</span>
<span id="L179" rel="#L179">179</span>
<span id="L180" rel="#L180">180</span>
<span id="L181" rel="#L181">181</span>
<span id="L182" rel="#L182">182</span>
<span id="L183" rel="#L183">183</span>
<span id="L184" rel="#L184">184</span>
<span id="L185" rel="#L185">185</span>
<span id="L186" rel="#L186">186</span>
<span id="L187" rel="#L187">187</span>
<span id="L188" rel="#L188">188</span>
<span id="L189" rel="#L189">189</span>
<span id="L190" rel="#L190">190</span>
<span id="L191" rel="#L191">191</span>
<span id="L192" rel="#L192">192</span>
<span id="L193" rel="#L193">193</span>
<span id="L194" rel="#L194">194</span>
<span id="L195" rel="#L195">195</span>
<span id="L196" rel="#L196">196</span>
<span id="L197" rel="#L197">197</span>
<span id="L198" rel="#L198">198</span>
<span id="L199" rel="#L199">199</span>
<span id="L200" rel="#L200">200</span>
<span id="L201" rel="#L201">201</span>
<span id="L202" rel="#L202">202</span>
<span id="L203" rel="#L203">203</span>
<span id="L204" rel="#L204">204</span>
<span id="L205" rel="#L205">205</span>
<span id="L206" rel="#L206">206</span>
<span id="L207" rel="#L207">207</span>
<span id="L208" rel="#L208">208</span>
<span id="L209" rel="#L209">209</span>
<span id="L210" rel="#L210">210</span>
<span id="L211" rel="#L211">211</span>
<span id="L212" rel="#L212">212</span>
<span id="L213" rel="#L213">213</span>
<span id="L214" rel="#L214">214</span>
<span id="L215" rel="#L215">215</span>
<span id="L216" rel="#L216">216</span>
<span id="L217" rel="#L217">217</span>
<span id="L218" rel="#L218">218</span>
<span id="L219" rel="#L219">219</span>
<span id="L220" rel="#L220">220</span>
<span id="L221" rel="#L221">221</span>
<span id="L222" rel="#L222">222</span>
<span id="L223" rel="#L223">223</span>
<span id="L224" rel="#L224">224</span>
<span id="L225" rel="#L225">225</span>
<span id="L226" rel="#L226">226</span>
<span id="L227" rel="#L227">227</span>
<span id="L228" rel="#L228">228</span>
<span id="L229" rel="#L229">229</span>
<span id="L230" rel="#L230">230</span>
<span id="L231" rel="#L231">231</span>
<span id="L232" rel="#L232">232</span>
<span id="L233" rel="#L233">233</span>
<span id="L234" rel="#L234">234</span>
<span id="L235" rel="#L235">235</span>
<span id="L236" rel="#L236">236</span>
<span id="L237" rel="#L237">237</span>
<span id="L238" rel="#L238">238</span>
<span id="L239" rel="#L239">239</span>
<span id="L240" rel="#L240">240</span>
<span id="L241" rel="#L241">241</span>
<span id="L242" rel="#L242">242</span>
<span id="L243" rel="#L243">243</span>
<span id="L244" rel="#L244">244</span>
<span id="L245" rel="#L245">245</span>
<span id="L246" rel="#L246">246</span>
<span id="L247" rel="#L247">247</span>
<span id="L248" rel="#L248">248</span>
<span id="L249" rel="#L249">249</span>
<span id="L250" rel="#L250">250</span>
<span id="L251" rel="#L251">251</span>
<span id="L252" rel="#L252">252</span>
<span id="L253" rel="#L253">253</span>
<span id="L254" rel="#L254">254</span>
<span id="L255" rel="#L255">255</span>
<span id="L256" rel="#L256">256</span>
<span id="L257" rel="#L257">257</span>
<span id="L258" rel="#L258">258</span>
<span id="L259" rel="#L259">259</span>
<span id="L260" rel="#L260">260</span>
<span id="L261" rel="#L261">261</span>
<span id="L262" rel="#L262">262</span>
<span id="L263" rel="#L263">263</span>
<span id="L264" rel="#L264">264</span>
<span id="L265" rel="#L265">265</span>
<span id="L266" rel="#L266">266</span>
<span id="L267" rel="#L267">267</span>
<span id="L268" rel="#L268">268</span>
<span id="L269" rel="#L269">269</span>
<span id="L270" rel="#L270">270</span>
<span id="L271" rel="#L271">271</span>
<span id="L272" rel="#L272">272</span>
<span id="L273" rel="#L273">273</span>
<span id="L274" rel="#L274">274</span>
<span id="L275" rel="#L275">275</span>
<span id="L276" rel="#L276">276</span>
<span id="L277" rel="#L277">277</span>
<span id="L278" rel="#L278">278</span>
<span id="L279" rel="#L279">279</span>
<span id="L280" rel="#L280">280</span>
<span id="L281" rel="#L281">281</span>
<span id="L282" rel="#L282">282</span>
<span id="L283" rel="#L283">283</span>
<span id="L284" rel="#L284">284</span>
<span id="L285" rel="#L285">285</span>
<span id="L286" rel="#L286">286</span>
<span id="L287" rel="#L287">287</span>
<span id="L288" rel="#L288">288</span>
<span id="L289" rel="#L289">289</span>
<span id="L290" rel="#L290">290</span>
<span id="L291" rel="#L291">291</span>
<span id="L292" rel="#L292">292</span>
<span id="L293" rel="#L293">293</span>
<span id="L294" rel="#L294">294</span>
<span id="L295" rel="#L295">295</span>
<span id="L296" rel="#L296">296</span>
<span id="L297" rel="#L297">297</span>
<span id="L298" rel="#L298">298</span>

            </td>
            <td class="blob-line-code">
                    <div class="highlight"><pre><div class='line' id='LC1'><span class="cm">/* The MIT License</span></div><div class='line' id='LC2'><br/></div><div class='line' id='LC3'><span class="cm">   Copyright (c) 2008, 2011 Attractive Chaos &lt;attractor@live.co.uk&gt;</span></div><div class='line' id='LC4'><br/></div><div class='line' id='LC5'><span class="cm">   Permission is hereby granted, free of charge, to any person obtaining</span></div><div class='line' id='LC6'><span class="cm">   a copy of this software and associated documentation files (the</span></div><div class='line' id='LC7'><span class="cm">   &quot;Software&quot;), to deal in the Software without restriction, including</span></div><div class='line' id='LC8'><span class="cm">   without limitation the rights to use, copy, modify, merge, publish,</span></div><div class='line' id='LC9'><span class="cm">   distribute, sublicense, and/or sell copies of the Software, and to</span></div><div class='line' id='LC10'><span class="cm">   permit persons to whom the Software is furnished to do so, subject to</span></div><div class='line' id='LC11'><span class="cm">   the following conditions:</span></div><div class='line' id='LC12'><br/></div><div class='line' id='LC13'><span class="cm">   The above copyright notice and this permission notice shall be</span></div><div class='line' id='LC14'><span class="cm">   included in all copies or substantial portions of the Software.</span></div><div class='line' id='LC15'><br/></div><div class='line' id='LC16'><span class="cm">   THE SOFTWARE IS PROVIDED &quot;AS IS&quot;, WITHOUT WARRANTY OF ANY KIND,</span></div><div class='line' id='LC17'><span class="cm">   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF</span></div><div class='line' id='LC18'><span class="cm">   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND</span></div><div class='line' id='LC19'><span class="cm">   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS</span></div><div class='line' id='LC20'><span class="cm">   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN</span></div><div class='line' id='LC21'><span class="cm">   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN</span></div><div class='line' id='LC22'><span class="cm">   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE</span></div><div class='line' id='LC23'><span class="cm">   SOFTWARE.</span></div><div class='line' id='LC24'><span class="cm">*/</span></div><div class='line' id='LC25'><br/></div><div class='line' id='LC26'><span class="cm">/*</span></div><div class='line' id='LC27'><span class="cm">  2011-04-10 (0.1.6):</span></div><div class='line' id='LC28'><br/></div><div class='line' id='LC29'><span class="cm">  	* Added sample</span></div><div class='line' id='LC30'><br/></div><div class='line' id='LC31'><span class="cm">  2011-03 (0.1.5):</span></div><div class='line' id='LC32'><br/></div><div class='line' id='LC33'><span class="cm">	* Added shuffle/permutation</span></div><div class='line' id='LC34'><br/></div><div class='line' id='LC35'><span class="cm">  2008-11-16 (0.1.4):</span></div><div class='line' id='LC36'><br/></div><div class='line' id='LC37'><span class="cm">    * Fixed a bug in introsort() that happens in rare cases.</span></div><div class='line' id='LC38'><br/></div><div class='line' id='LC39'><span class="cm">  2008-11-05 (0.1.3):</span></div><div class='line' id='LC40'><br/></div><div class='line' id='LC41'><span class="cm">    * Fixed a bug in introsort() for complex comparisons.</span></div><div class='line' id='LC42'><br/></div><div class='line' id='LC43'><span class="cm">	* Fixed a bug in mergesort(). The previous version is not stable.</span></div><div class='line' id='LC44'><br/></div><div class='line' id='LC45'><span class="cm">  2008-09-15 (0.1.2):</span></div><div class='line' id='LC46'><br/></div><div class='line' id='LC47'><span class="cm">	* Accelerated introsort. On my Mac (not on another Linux machine),</span></div><div class='line' id='LC48'><span class="cm">	  my implementation is as fast as std::sort on random input.</span></div><div class='line' id='LC49'><br/></div><div class='line' id='LC50'><span class="cm">	* Added combsort and in introsort, switch to combsort if the</span></div><div class='line' id='LC51'><span class="cm">	  recursion is too deep.</span></div><div class='line' id='LC52'><br/></div><div class='line' id='LC53'><span class="cm">  2008-09-13 (0.1.1):</span></div><div class='line' id='LC54'><br/></div><div class='line' id='LC55'><span class="cm">	* Added k-small algorithm</span></div><div class='line' id='LC56'><br/></div><div class='line' id='LC57'><span class="cm">  2008-09-05 (0.1.0):</span></div><div class='line' id='LC58'><br/></div><div class='line' id='LC59'><span class="cm">	* Initial version</span></div><div class='line' id='LC60'><br/></div><div class='line' id='LC61'><span class="cm">*/</span></div><div class='line' id='LC62'><br/></div><div class='line' id='LC63'><span class="cp">#ifndef AC_KSORT_H</span></div><div class='line' id='LC64'><span class="cp">#define AC_KSORT_H</span></div><div class='line' id='LC65'><br/></div><div class='line' id='LC66'><span class="cp">#include &lt;stdlib.h&gt;</span></div><div class='line' id='LC67'><span class="cp">#include &lt;string.h&gt;</span></div><div class='line' id='LC68'><br/></div><div class='line' id='LC69'><span class="k">typedef</span> <span class="k">struct</span> <span class="p">{</span></div><div class='line' id='LC70'>	<span class="kt">void</span> <span class="o">*</span><span class="n">left</span><span class="p">,</span> <span class="o">*</span><span class="n">right</span><span class="p">;</span></div><div class='line' id='LC71'>	<span class="kt">int</span> <span class="n">depth</span><span class="p">;</span></div><div class='line' id='LC72'><span class="p">}</span> <span class="kt">ks_isort_stack_t</span><span class="p">;</span></div><div class='line' id='LC73'><br/></div><div class='line' id='LC74'><span class="cp">#define KSORT_SWAP(type_t, a, b) { register type_t t=(a); (a)=(b); (b)=t; }</span></div><div class='line' id='LC75'><br/></div><div class='line' id='LC76'><span class="cp">#define KSORT_INIT(name, type_t, __sort_lt)								\</span></div><div class='line' id='LC77'><span class="cp">	void ks_mergesort_##name(size_t n, type_t array[], type_t temp[])	\</span></div><div class='line' id='LC78'><span class="cp">	{																	\</span></div><div class='line' id='LC79'><span class="cp">		type_t *a2[2], *a, *b;											\</span></div><div class='line' id='LC80'><span class="cp">		int curr, shift;												\</span></div><div class='line' id='LC81'><span class="cp">																		\</span></div><div class='line' id='LC82'><span class="cp">		a2[0] = array;													\</span></div><div class='line' id='LC83'><span class="cp">		a2[1] = temp? temp : (type_t*)malloc(sizeof(type_t) * n);		\</span></div><div class='line' id='LC84'><span class="cp">		for (curr = 0, shift = 0; (1ul&lt;&lt;shift) &lt; n; ++shift) {			\</span></div><div class='line' id='LC85'><span class="cp">			a = a2[curr]; b = a2[1-curr];								\</span></div><div class='line' id='LC86'><span class="cp">			if (shift == 0) {											\</span></div><div class='line' id='LC87'><span class="cp">				type_t *p = b, *i, *eb = a + n;							\</span></div><div class='line' id='LC88'><span class="cp">				for (i = a; i &lt; eb; i += 2) {							\</span></div><div class='line' id='LC89'><span class="cp">					if (i == eb - 1) *p++ = *i;							\</span></div><div class='line' id='LC90'><span class="cp">					else {												\</span></div><div class='line' id='LC91'><span class="cp">						if (__sort_lt(*(i+1), *i)) {					\</span></div><div class='line' id='LC92'><span class="cp">							*p++ = *(i+1); *p++ = *i;					\</span></div><div class='line' id='LC93'><span class="cp">						} else {										\</span></div><div class='line' id='LC94'><span class="cp">							*p++ = *i; *p++ = *(i+1);					\</span></div><div class='line' id='LC95'><span class="cp">						}												\</span></div><div class='line' id='LC96'><span class="cp">					}													\</span></div><div class='line' id='LC97'><span class="cp">				}														\</span></div><div class='line' id='LC98'><span class="cp">			} else {													\</span></div><div class='line' id='LC99'><span class="cp">				size_t i, step = 1ul&lt;&lt;shift;							\</span></div><div class='line' id='LC100'><span class="cp">				for (i = 0; i &lt; n; i += step&lt;&lt;1) {						\</span></div><div class='line' id='LC101'><span class="cp">					type_t *p, *j, *k, *ea, *eb;						\</span></div><div class='line' id='LC102'><span class="cp">					if (n &lt; i + step) {									\</span></div><div class='line' id='LC103'><span class="cp">						ea = a + n; eb = a;								\</span></div><div class='line' id='LC104'><span class="cp">					} else {											\</span></div><div class='line' id='LC105'><span class="cp">						ea = a + i + step;								\</span></div><div class='line' id='LC106'><span class="cp">						eb = a + (n &lt; i + (step&lt;&lt;1)? n : i + (step&lt;&lt;1)); \</span></div><div class='line' id='LC107'><span class="cp">					}													\</span></div><div class='line' id='LC108'><span class="cp">					j = a + i; k = a + i + step; p = b + i;				\</span></div><div class='line' id='LC109'><span class="cp">					while (j &lt; ea &amp;&amp; k &lt; eb) {							\</span></div><div class='line' id='LC110'><span class="cp">						if (__sort_lt(*k, *j)) *p++ = *k++;				\</span></div><div class='line' id='LC111'><span class="cp">						else *p++ = *j++;								\</span></div><div class='line' id='LC112'><span class="cp">					}													\</span></div><div class='line' id='LC113'><span class="cp">					while (j &lt; ea) *p++ = *j++;							\</span></div><div class='line' id='LC114'><span class="cp">					while (k &lt; eb) *p++ = *k++;							\</span></div><div class='line' id='LC115'><span class="cp">				}														\</span></div><div class='line' id='LC116'><span class="cp">			}															\</span></div><div class='line' id='LC117'><span class="cp">			curr = 1 - curr;											\</span></div><div class='line' id='LC118'><span class="cp">		}																\</span></div><div class='line' id='LC119'><span class="cp">		if (curr == 1) {												\</span></div><div class='line' id='LC120'><span class="cp">			type_t *p = a2[0], *i = a2[1], *eb = array + n;				\</span></div><div class='line' id='LC121'><span class="cp">			for (; p &lt; eb; ++i) *p++ = *i;								\</span></div><div class='line' id='LC122'><span class="cp">		}																\</span></div><div class='line' id='LC123'><span class="cp">		if (temp == 0) free(a2[1]);										\</span></div><div class='line' id='LC124'><span class="cp">	}																	\</span></div><div class='line' id='LC125'><span class="cp">	void ks_heapadjust_##name(size_t i, size_t n, type_t l[])			\</span></div><div class='line' id='LC126'><span class="cp">	{																	\</span></div><div class='line' id='LC127'><span class="cp">		size_t k = i;													\</span></div><div class='line' id='LC128'><span class="cp">		type_t tmp = l[i];												\</span></div><div class='line' id='LC129'><span class="cp">		while ((k = (k &lt;&lt; 1) + 1) &lt; n) {								\</span></div><div class='line' id='LC130'><span class="cp">			if (k != n - 1 &amp;&amp; __sort_lt(l[k], l[k+1])) ++k;				\</span></div><div class='line' id='LC131'><span class="cp">			if (__sort_lt(l[k], tmp)) break;							\</span></div><div class='line' id='LC132'><span class="cp">			l[i] = l[k]; i = k;											\</span></div><div class='line' id='LC133'><span class="cp">		}																\</span></div><div class='line' id='LC134'><span class="cp">		l[i] = tmp;														\</span></div><div class='line' id='LC135'><span class="cp">	}																	\</span></div><div class='line' id='LC136'><span class="cp">	void ks_heapmake_##name(size_t lsize, type_t l[])					\</span></div><div class='line' id='LC137'><span class="cp">	{																	\</span></div><div class='line' id='LC138'><span class="cp">		size_t i;														\</span></div><div class='line' id='LC139'><span class="cp">		for (i = (lsize &gt;&gt; 1) - 1; i != (size_t)(-1); --i)				\</span></div><div class='line' id='LC140'><span class="cp">			ks_heapadjust_##name(i, lsize, l);							\</span></div><div class='line' id='LC141'><span class="cp">	}																	\</span></div><div class='line' id='LC142'><span class="cp">	void ks_heapsort_##name(size_t lsize, type_t l[])					\</span></div><div class='line' id='LC143'><span class="cp">	{																	\</span></div><div class='line' id='LC144'><span class="cp">		size_t i;														\</span></div><div class='line' id='LC145'><span class="cp">		for (i = lsize - 1; i &gt; 0; --i) {								\</span></div><div class='line' id='LC146'><span class="cp">			type_t tmp;													\</span></div><div class='line' id='LC147'><span class="cp">			tmp = *l; *l = l[i]; l[i] = tmp; ks_heapadjust_##name(0, i, l); \</span></div><div class='line' id='LC148'><span class="cp">		}																\</span></div><div class='line' id='LC149'><span class="cp">	}																	\</span></div><div class='line' id='LC150'><span class="cp">	inline void __ks_insertsort_##name(type_t *s, type_t *t)			\</span></div><div class='line' id='LC151'><span class="cp">	{																	\</span></div><div class='line' id='LC152'><span class="cp">		type_t *i, *j, swap_tmp;										\</span></div><div class='line' id='LC153'><span class="cp">		for (i = s + 1; i &lt; t; ++i)										\</span></div><div class='line' id='LC154'><span class="cp">			for (j = i; j &gt; s &amp;&amp; __sort_lt(*j, *(j-1)); --j) {			\</span></div><div class='line' id='LC155'><span class="cp">				swap_tmp = *j; *j = *(j-1); *(j-1) = swap_tmp;			\</span></div><div class='line' id='LC156'><span class="cp">			}															\</span></div><div class='line' id='LC157'><span class="cp">	}																	\</span></div><div class='line' id='LC158'><span class="cp">	void ks_combsort_##name(size_t n, type_t a[])						\</span></div><div class='line' id='LC159'><span class="cp">	{																	\</span></div><div class='line' id='LC160'><span class="cp">		const double shrink_factor = 1.2473309501039786540366528676643; \</span></div><div class='line' id='LC161'><span class="cp">		int do_swap;													\</span></div><div class='line' id='LC162'><span class="cp">		size_t gap = n;													\</span></div><div class='line' id='LC163'><span class="cp">		type_t tmp, *i, *j;												\</span></div><div class='line' id='LC164'><span class="cp">		do {															\</span></div><div class='line' id='LC165'><span class="cp">			if (gap &gt; 2) {												\</span></div><div class='line' id='LC166'><span class="cp">				gap = (size_t)(gap / shrink_factor);					\</span></div><div class='line' id='LC167'><span class="cp">				if (gap == 9 || gap == 10) gap = 11;					\</span></div><div class='line' id='LC168'><span class="cp">			}															\</span></div><div class='line' id='LC169'><span class="cp">			do_swap = 0;												\</span></div><div class='line' id='LC170'><span class="cp">			for (i = a; i &lt; a + n - gap; ++i) {							\</span></div><div class='line' id='LC171'><span class="cp">				j = i + gap;											\</span></div><div class='line' id='LC172'><span class="cp">				if (__sort_lt(*j, *i)) {								\</span></div><div class='line' id='LC173'><span class="cp">					tmp = *i; *i = *j; *j = tmp;						\</span></div><div class='line' id='LC174'><span class="cp">					do_swap = 1;										\</span></div><div class='line' id='LC175'><span class="cp">				}														\</span></div><div class='line' id='LC176'><span class="cp">			}															\</span></div><div class='line' id='LC177'><span class="cp">		} while (do_swap || gap &gt; 2);									\</span></div><div class='line' id='LC178'><span class="cp">		if (gap != 1) __ks_insertsort_##name(a, a + n);					\</span></div><div class='line' id='LC179'><span class="cp">	}																	\</span></div><div class='line' id='LC180'><span class="cp">	void ks_introsort_##name(size_t n, type_t a[])						\</span></div><div class='line' id='LC181'><span class="cp">	{																	\</span></div><div class='line' id='LC182'><span class="cp">		int d;															\</span></div><div class='line' id='LC183'><span class="cp">		ks_isort_stack_t *top, *stack;									\</span></div><div class='line' id='LC184'><span class="cp">		type_t rp, swap_tmp;											\</span></div><div class='line' id='LC185'><span class="cp">		type_t *s, *t, *i, *j, *k;										\</span></div><div class='line' id='LC186'><span class="cp">																		\</span></div><div class='line' id='LC187'><span class="cp">		if (n &lt; 1) return;												\</span></div><div class='line' id='LC188'><span class="cp">		else if (n == 2) {												\</span></div><div class='line' id='LC189'><span class="cp">			if (__sort_lt(a[1], a[0])) { swap_tmp = a[0]; a[0] = a[1]; a[1] = swap_tmp; } \</span></div><div class='line' id='LC190'><span class="cp">			return;														\</span></div><div class='line' id='LC191'><span class="cp">		}																\</span></div><div class='line' id='LC192'><span class="cp">		for (d = 2; 1ul&lt;&lt;d &lt; n; ++d);									\</span></div><div class='line' id='LC193'><span class="cp">		stack = (ks_isort_stack_t*)malloc(sizeof(ks_isort_stack_t) * ((sizeof(size_t)*d)+2)); \</span></div><div class='line' id='LC194'><span class="cp">		top = stack; s = a; t = a + (n-1); d &lt;&lt;= 1;						\</span></div><div class='line' id='LC195'><span class="cp">		while (1) {														\</span></div><div class='line' id='LC196'><span class="cp">			if (s &lt; t) {												\</span></div><div class='line' id='LC197'><span class="cp">				if (--d == 0) {											\</span></div><div class='line' id='LC198'><span class="cp">					ks_combsort_##name(t - s + 1, s);					\</span></div><div class='line' id='LC199'><span class="cp">					t = s;												\</span></div><div class='line' id='LC200'><span class="cp">					continue;											\</span></div><div class='line' id='LC201'><span class="cp">				}														\</span></div><div class='line' id='LC202'><span class="cp">				i = s; j = t; k = i + ((j-i)&gt;&gt;1) + 1;					\</span></div><div class='line' id='LC203'><span class="cp">				if (__sort_lt(*k, *i)) {								\</span></div><div class='line' id='LC204'><span class="cp">					if (__sort_lt(*k, *j)) k = j;						\</span></div><div class='line' id='LC205'><span class="cp">				} else k = __sort_lt(*j, *i)? i : j;					\</span></div><div class='line' id='LC206'><span class="cp">				rp = *k;												\</span></div><div class='line' id='LC207'><span class="cp">				if (k != t) { swap_tmp = *k; *k = *t; *t = swap_tmp; }	\</span></div><div class='line' id='LC208'><span class="cp">				for (;;) {												\</span></div><div class='line' id='LC209'><span class="cp">					do ++i; while (__sort_lt(*i, rp));					\</span></div><div class='line' id='LC210'><span class="cp">					do --j; while (i &lt;= j &amp;&amp; __sort_lt(rp, *j));		\</span></div><div class='line' id='LC211'><span class="cp">					if (j &lt;= i) break;									\</span></div><div class='line' id='LC212'><span class="cp">					swap_tmp = *i; *i = *j; *j = swap_tmp;				\</span></div><div class='line' id='LC213'><span class="cp">				}														\</span></div><div class='line' id='LC214'><span class="cp">				swap_tmp = *i; *i = *t; *t = swap_tmp;					\</span></div><div class='line' id='LC215'><span class="cp">				if (i-s &gt; t-i) {										\</span></div><div class='line' id='LC216'><span class="cp">					if (i-s &gt; 16) { top-&gt;left = s; top-&gt;right = i-1; top-&gt;depth = d; ++top; } \</span></div><div class='line' id='LC217'><span class="cp">					s = t-i &gt; 16? i+1 : t;								\</span></div><div class='line' id='LC218'><span class="cp">				} else {												\</span></div><div class='line' id='LC219'><span class="cp">					if (t-i &gt; 16) { top-&gt;left = i+1; top-&gt;right = t; top-&gt;depth = d; ++top; } \</span></div><div class='line' id='LC220'><span class="cp">					t = i-s &gt; 16? i-1 : s;								\</span></div><div class='line' id='LC221'><span class="cp">				}														\</span></div><div class='line' id='LC222'><span class="cp">			} else {													\</span></div><div class='line' id='LC223'><span class="cp">				if (top == stack) {										\</span></div><div class='line' id='LC224'><span class="cp">					free(stack);										\</span></div><div class='line' id='LC225'><span class="cp">					__ks_insertsort_##name(a, a+n);						\</span></div><div class='line' id='LC226'><span class="cp">					return;												\</span></div><div class='line' id='LC227'><span class="cp">				} else { --top; s = (type_t*)top-&gt;left; t = (type_t*)top-&gt;right; d = top-&gt;depth; } \</span></div><div class='line' id='LC228'><span class="cp">			}															\</span></div><div class='line' id='LC229'><span class="cp">		}																\</span></div><div class='line' id='LC230'><span class="cp">	}																	\</span></div><div class='line' id='LC231'><span class="cp">	</span><span class="cm">/* This function is adapted from: http://ndevilla.free.fr/median/ */</span><span class="cp"> \</span></div><div class='line' id='LC232'><span class="cp">	</span><span class="cm">/* 0 &lt;= kk &lt; n */</span><span class="cp">													\</span></div><div class='line' id='LC233'><span class="cp">	type_t ks_ksmall_##name(size_t n, type_t arr[], size_t kk)			\</span></div><div class='line' id='LC234'><span class="cp">	{																	\</span></div><div class='line' id='LC235'><span class="cp">		type_t *low, *high, *k, *ll, *hh, *mid;							\</span></div><div class='line' id='LC236'><span class="cp">		low = arr; high = arr + n - 1; k = arr + kk;					\</span></div><div class='line' id='LC237'><span class="cp">		for (;;) {														\</span></div><div class='line' id='LC238'><span class="cp">			if (high &lt;= low) return *k;									\</span></div><div class='line' id='LC239'><span class="cp">			if (high == low + 1) {										\</span></div><div class='line' id='LC240'><span class="cp">				if (__sort_lt(*high, *low)) KSORT_SWAP(type_t, *low, *high); \</span></div><div class='line' id='LC241'><span class="cp">				return *k;												\</span></div><div class='line' id='LC242'><span class="cp">			}															\</span></div><div class='line' id='LC243'><span class="cp">			mid = low + (high - low) / 2;								\</span></div><div class='line' id='LC244'><span class="cp">			if (__sort_lt(*high, *mid)) KSORT_SWAP(type_t, *mid, *high); \</span></div><div class='line' id='LC245'><span class="cp">			if (__sort_lt(*high, *low)) KSORT_SWAP(type_t, *low, *high); \</span></div><div class='line' id='LC246'><span class="cp">			if (__sort_lt(*low, *mid)) KSORT_SWAP(type_t, *mid, *low);	\</span></div><div class='line' id='LC247'><span class="cp">			KSORT_SWAP(type_t, *mid, *(low+1));							\</span></div><div class='line' id='LC248'><span class="cp">			ll = low + 1; hh = high;									\</span></div><div class='line' id='LC249'><span class="cp">			for (;;) {													\</span></div><div class='line' id='LC250'><span class="cp">				do ++ll; while (__sort_lt(*ll, *low));					\</span></div><div class='line' id='LC251'><span class="cp">				do --hh; while (__sort_lt(*low, *hh));					\</span></div><div class='line' id='LC252'><span class="cp">				if (hh &lt; ll) break;										\</span></div><div class='line' id='LC253'><span class="cp">				KSORT_SWAP(type_t, *ll, *hh);							\</span></div><div class='line' id='LC254'><span class="cp">			}															\</span></div><div class='line' id='LC255'><span class="cp">			KSORT_SWAP(type_t, *low, *hh);								\</span></div><div class='line' id='LC256'><span class="cp">			if (hh &lt;= k) low = ll;										\</span></div><div class='line' id='LC257'><span class="cp">			if (hh &gt;= k) high = hh - 1;									\</span></div><div class='line' id='LC258'><span class="cp">		}																\</span></div><div class='line' id='LC259'><span class="cp">	}																	\</span></div><div class='line' id='LC260'><span class="cp">	void ks_shuffle_##name(size_t n, type_t a[])						\</span></div><div class='line' id='LC261'><span class="cp">	{																	\</span></div><div class='line' id='LC262'><span class="cp">		int i, j;														\</span></div><div class='line' id='LC263'><span class="cp">		for (i = n; i &gt; 1; --i) {										\</span></div><div class='line' id='LC264'><span class="cp">			type_t tmp;													\</span></div><div class='line' id='LC265'><span class="cp">			j = (int)(drand48() * i);									\</span></div><div class='line' id='LC266'><span class="cp">			tmp = a[j]; a[j] = a[i-1]; a[i-1] = tmp;					\</span></div><div class='line' id='LC267'><span class="cp">		}																\</span></div><div class='line' id='LC268'><span class="cp">	}																	\</span></div><div class='line' id='LC269'><span class="cp">	void ks_sample_##name(size_t n, size_t r, type_t a[]) </span><span class="cm">/* FIXME: NOT TESTED!!! */</span><span class="cp"> \</span></div><div class='line' id='LC270'><span class="cp">	{ </span><span class="cm">/* reference: http://code.activestate.com/recipes/272884/ */</span><span class="cp"> \</span></div><div class='line' id='LC271'><span class="cp">		int i, k, pop = n; \</span></div><div class='line' id='LC272'><span class="cp">		for (i = (int)r, k = 0; i &gt;= 0; --i) { \</span></div><div class='line' id='LC273'><span class="cp">			double z = 1., x = drand48(); \</span></div><div class='line' id='LC274'><span class="cp">			type_t tmp; \</span></div><div class='line' id='LC275'><span class="cp">			while (x &lt; z) z -= z * i / (pop--); \</span></div><div class='line' id='LC276'><span class="cp">			if (k != n - pop - 1) tmp = a[k], a[k] = a[n-pop-1], a[n-pop-1] = tmp; \</span></div><div class='line' id='LC277'><span class="cp">			++k; \</span></div><div class='line' id='LC278'><span class="cp">		} \</span></div><div class='line' id='LC279'><span class="cp">	}</span></div><div class='line' id='LC280'><br/></div><div class='line' id='LC281'><span class="cp">#define ks_mergesort(name, n, a, t) ks_mergesort_##name(n, a, t)</span></div><div class='line' id='LC282'><span class="cp">#define ks_introsort(name, n, a) ks_introsort_##name(n, a)</span></div><div class='line' id='LC283'><span class="cp">#define ks_combsort(name, n, a) ks_combsort_##name(n, a)</span></div><div class='line' id='LC284'><span class="cp">#define ks_heapsort(name, n, a) ks_heapsort_##name(n, a)</span></div><div class='line' id='LC285'><span class="cp">#define ks_heapmake(name, n, a) ks_heapmake_##name(n, a)</span></div><div class='line' id='LC286'><span class="cp">#define ks_heapadjust(name, i, n, a) ks_heapadjust_##name(i, n, a)</span></div><div class='line' id='LC287'><span class="cp">#define ks_ksmall(name, n, a, k) ks_ksmall_##name(n, a, k)</span></div><div class='line' id='LC288'><span class="cp">#define ks_shuffle(name, n, a) ks_shuffle_##name(n, a)</span></div><div class='line' id='LC289'><br/></div><div class='line' id='LC290'><span class="cp">#define ks_lt_generic(a, b) ((a) &lt; (b))</span></div><div class='line' id='LC291'><span class="cp">#define ks_lt_str(a, b) (strcmp((a), (b)) &lt; 0)</span></div><div class='line' id='LC292'><br/></div><div class='line' id='LC293'><span class="k">typedef</span> <span class="k">const</span> <span class="kt">char</span> <span class="o">*</span><span class="kt">ksstr_t</span><span class="p">;</span></div><div class='line' id='LC294'><br/></div><div class='line' id='LC295'><span class="cp">#define KSORT_INIT_GENERIC(type_t) KSORT_INIT(type_t, type_t, ks_lt_generic)</span></div><div class='line' id='LC296'><span class="cp">#define KSORT_INIT_STR KSORT_INIT(str, ksstr_t, ks_lt_str)</span></div><div class='line' id='LC297'><br/></div><div class='line' id='LC298'><span class="cp">#endif</span></div></pre></div>
            </td>
          </tr>
        </table>
  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2013 <span title="0.06067s from fe2.rs.github.com">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
          <div class="suggester-container">
              <div class="suggester fullscreen-suggester js-navigation-container" id="fullscreen_suggester"
                 data-url="/attractivechaos/klib/suggestions/commit">
              </div>
          </div>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped leftwards" title="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped leftwards"
      title="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-remove-close close ajax-error-dismiss"></a>
      Something went wrong with that request. Please try again.
    </div>

    
  </body>
</html>

