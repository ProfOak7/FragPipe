package com.dmtavt.fragpipe.tabs;

import com.dmtavt.fragpipe.params.ptmprophet.PtmProphetPanel;
import com.github.chhh.utils.swing.JPanelWithEnablement;
import com.github.chhh.utils.swing.MigUtils;
import com.dmtavt.fragpipe.tools.ptmshepherd.PtmshepherdPanel;

public class TabPtms extends JPanelWithEnablement {
  private static MigUtils mu = MigUtils.get();
  public static final String TAB_PREFIX = "quantitation.";
  private PtmshepherdPanel panelPtmshepherd;
  private PtmProphetPanel panelPtmProhpet;

  public TabPtms() {
    init();
    initMore();
  }

  private void initMore() {
    //Bus.registerQuietly(this);
  }

  private void init() {
    mu.layout(this).fillX();

    panelPtmshepherd = new PtmshepherdPanel();

    mu.add(this, panelPtmshepherd).growX().wrap();
  }
}
