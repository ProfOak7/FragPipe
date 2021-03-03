package com.dmtavt.fragpipe.tabs;

import com.dmtavt.fragpipe.tools.pepproph.PepProphPanel;
import com.dmtavt.fragpipe.tools.percolator.PercolatorPanel;
import com.dmtavt.fragpipe.tools.protproph.ProtProphPanel;
import com.dmtavt.fragpipe.tools.ptmprophet.PtmProphetPanel;
import com.github.chhh.utils.swing.JPanelWithEnablement;
import com.github.chhh.utils.swing.MigUtils;
import com.dmtavt.fragpipe.tools.crystalc.CrystalcPanel;
import com.dmtavt.fragpipe.tools.philosopher.ReportPanel;

public class TabValidation extends JPanelWithEnablement {
  private static MigUtils mu = MigUtils.get();
  public static final String TAB_PREFIX = "validation.";
  private PercolatorPanel panelPercolator;
  private PepProphPanel panelPepProph;
  private ProtProphPanel panelProtProph;
  private CrystalcPanel panelCrystalc;
  private ReportPanel panelReport;
  private PtmProphetPanel panelPtmProphet;

  public TabValidation() {
    init();
    initMore();
  }

  private void initMore() {
    //Bus.registerQuietly(this);
  }

  private void init() {
    mu.layout(this).fillX();

    panelCrystalc = new CrystalcPanel();
    panelPercolator = new PercolatorPanel();
    panelPepProph = new PepProphPanel();
    panelPtmProphet = new PtmProphetPanel();
    panelProtProph = new ProtProphPanel();
    panelReport = new ReportPanel();


    mu.add(this, panelCrystalc).growX().wrap();
    mu.add(this, panelPercolator).growX().wrap();
    mu.add(this, panelPepProph).growX().wrap();
    mu.add(this, panelPtmProphet).growX().wrap();
    mu.add(this, panelProtProph).growX().wrap();
    mu.add(this, panelReport).growX().wrap();
  }


}
