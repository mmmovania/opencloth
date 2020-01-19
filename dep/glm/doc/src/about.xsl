<?xml version="1.0" encoding="iso-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="html" media-type="text/html; charset=ISO-8859-1" />

  <xsl:include href="./util.xsl" />

  <xsl:template match="/">
    <html>
      <head>
        <title>OpenGL Mathematics : About</title>
        <meta http-equiv="Content-Language" content="en" />
        <meta http-equiv="Content-Type" content="application/xhtml+xml; charset=iso-8859-1" />
        <meta name="copyright" content="G-Truc Creation" />
        <link href="./common/style.css" rel="stylesheet" media="screen, print, embossed" type="text/css" />
      </head>
      <body>
        <table>
          <tr>
            <td class="menu">
              <xsl:apply-templates select="./glm/menu" />
            </td>
            <td class="page">
              <div class="title1">
                <img src="./common/title.png" alt="OpenGL Mathematics" />
              </div>
              <xsl:apply-templates select="./glm/about-short" />
              <br />
              <xsl:apply-templates select="./glm/about-long" />

              <div class="email">
                <img src="./common/email.png" alt="email not available as text" />
              </div>
              <div class="news-separator">_________________</div>
              <br />
              <div class="title3">
                <xsl:value-of select="./glm/@copyright" />
                <a href="http://www.g-truc.net">G-Truc Creation</a>
              </div>
            </td>
          </tr>
        </table>
      </body>
    </html>
  </xsl:template>

  <xsl:template match="about-long">
    <div>
      <div class="title-date">
        <xsl:value-of select="./@date" />
      </div>
      <div class="title4">
        <xsl:value-of select="./@title" />
      </div>
      <div>
        <xsl:if test="./paragraph">
          <xsl:apply-templates select="./paragraph" />
        </xsl:if>
        <xsl:if test="./list">
          <xsl:apply-templates select="./list" />
        </xsl:if>
        <xsl:apply-templates select="./source" />
      </div>
      <div class="news-separator">_________________</div>
      <br />
    </div>
  </xsl:template>

</xsl:stylesheet>
