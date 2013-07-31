<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output method="text"/>

  <xsl:template match="/">
    <xsl:apply-templates select="region/germline" />
  </xsl:template>

  <xsl:template match="germline">
&gt;<xsl:value-of select="@name" /> <xsl:value-of select="@familyname" /> <xsl:value-of select="@cdr3index" />
<xsl:value-of select="@sequence" />
  </xsl:template>

</xsl:stylesheet>
