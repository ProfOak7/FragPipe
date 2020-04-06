/*
 * Copyright (C) 2018 Dmitry Avtonomov
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package com.github.chhh.utils;

import com.github.chhh.utils.swing.GhostedTextComponent;
import com.github.chhh.utils.swing.StringRepresentable;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Desktop;
import java.awt.Dialog;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.Window;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.function.Consumer;
import java.util.function.Predicate;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JToggleButton;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.HyperlinkEvent;
import javax.swing.table.DefaultTableModel;
import javax.swing.text.Document;
import javax.swing.text.JTextComponent;
import net.miginfocom.layout.LC;
import net.miginfocom.swing.MigLayout;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import umich.msfragger.gui.MsfraggerGuiFrame;

/**
 * @author dmitriya
 */
public class SwingUtils {
  private static final Logger log = LoggerFactory.getLogger(SwingUtils.class);
  private static volatile String[] fontNames = null;
  private static volatile Font[] fonts = null;
  private static final Object fontLock = new Object();

  private SwingUtils() {
  }

  /**
   * Wraps a component in JScrollPane and sets scroll-bar speed to a reasonable value.
   */
  public static JScrollPane scroll(Component comp) {
    JScrollPane s = new JScrollPane(comp);
    s.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
    s.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
    s.getVerticalScrollBar().setUnitIncrement(16);
    return s;
  }

  public static DialogAndThread runThreadWithProgressBar(String title, Component parent, Runnable runnable) {
    JFrame frame = SwingUtils.findParentFrame(parent);
    final JDialog dialog = new JDialog(frame, title, true);
    JProgressBar bar = new JProgressBar(0, 100);
    bar.setIndeterminate(true);
    Dimension d = new Dimension(300, 75);
    bar.setMinimumSize(d);
    bar.setSize(d);
    dialog.add(bar, BorderLayout.CENTER);
    dialog.setSize(d);
    dialog.setLocationRelativeTo(parent);

    Thread thread = new Thread(() -> {
      try {
        runnable.run();
      } catch (Exception ex) {
        throw new IllegalStateException("Something happened while running behind a progress bar", ex);
      } finally {
        dialog.setVisible(false);
        dialog.dispose();
      }

    });
    return new DialogAndThread(dialog, thread);
  }

  public static void setLaf() {
    /* Set the Nimbus look and feel */
    //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
    /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
     * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html
     */
    try {
      if (OsUtils.isWindows()) {
        // native look on windows
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
      } else {
        // nimbus otherwise
        for (UIManager.LookAndFeelInfo info : UIManager
            .getInstalledLookAndFeels()) {
          if ("Nimbus".equals(info.getName())) {
            UIManager.setLookAndFeel(info.getClassName());
            break;
          }
        }
      }
    } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException e1) {
      java.util.logging.Logger.getLogger(MsfraggerGuiFrame.class.getName())
          .log(java.util.logging.Level.SEVERE, null, e1);
      try {
        for (UIManager.LookAndFeelInfo info : UIManager
            .getInstalledLookAndFeels()) {
          if ("Nimbus".equals(info.getName())) {
            UIManager.setLookAndFeel(info.getClassName());
            break;
          }
        }
      } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException e2) {
        java.util.logging.Logger.getLogger(MsfraggerGuiFrame.class.getName())
            .log(java.util.logging.Level.SEVERE, null, e2);
      }
    }
    //</editor-fold>
  }

  public static class DialogAndThread {
    public final JDialog dialog;
    public final Thread thread;

    public DialogAndThread(JDialog dialog, Thread thread) {
      this.dialog = dialog;
      this.thread = thread;
    }
  }

  public static JTable tableFromTwoSiblingFiles(Map<Path, Path> paths) {
    String[] columns = {"From", "To", "At"};
    String[][] data = new String[paths.size()][3];
    int index = -1;
    for (Entry<Path, Path> kv : paths.entrySet()) {
      data[++index][0] = kv.getKey().getFileName().toString();
      data[index][1] = kv.getValue().getFileName().toString();
      if (!kv.getValue().getParent().equals(kv.getKey().getParent())) {
        throw new IllegalArgumentException("Files must be siblings");
      }
      data[index][2] = kv.getKey().getParent().toString();
    }

    DefaultTableModel model = new DefaultTableModel(data, columns);
    JTable table = new JTable(model);
    table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);

    return table;
  }

  public static String getStrVal(Component c) {

    String val;
    if (c instanceof StringRepresentable) {
      val = ((StringRepresentable) c).asString();
    } else if (c instanceof JFormattedTextField) {
      val = ((JFormattedTextField) c).getText();
    } else if (c instanceof JTextField) {
      val = ((JTextField) c).getText();
    } else if (c instanceof JSpinner) {
      val = ((JSpinner) c).getValue().toString();
    } else if (c instanceof JCheckBox) {
      val = Boolean.valueOf(((JCheckBox) c).isSelected()).toString();
    } else if (c instanceof JComboBox) {
      val = ((JComboBox<?>) c).getModel().getSelectedItem().toString();
    } else {
      throw new UnsupportedOperationException("getStrVal() not implemented for type: " + c.getClass().getCanonicalName());
    }

    return val.trim();
  }

  public static void setStrVal(Component c, String val) {
    if (c instanceof StringRepresentable) {
      ((StringRepresentable) c).fromString(val);
    } else if (c instanceof JFormattedTextField) {
      ((JFormattedTextField) c).setText(val);
    } else if (c instanceof JTextField) {
      ((JTextField) c).setText(val);
    } else if (c instanceof JCheckBox) {
      ((JCheckBox) c).setSelected(Boolean.parseBoolean(val));
    } else if (c instanceof JComboBox) {
      ((JComboBox<?>) c).getModel().setSelectedItem(val);
    } else if (c instanceof JSpinner) {
      ((JSpinner) c).setValue(Double.parseDouble(val));
    } else {
      throw new UnsupportedOperationException("setStrVal() not implemented for type: " + c.getClass().getCanonicalName());
    }
  }

  /**
   * Installs a listener to receive notification when the text of any {@code JTextComponent} is
   * changed. Internally, it installs a {@link DocumentListener} on the text component's {@link
   * Document}, and a {@link PropertyChangeListener} on the text component to detect if the {@code
   * Document} itself is replaced.
   *
   * @param text any text component, such as a {@link JTextField} or {@link JTextArea}
   * @param changeListener a listener to receieve {@link ChangeEvent}s when the text is changed; the
   * source object for the events will be the text component
   * @throws NullPointerException if either parameter is null
   *
   * Taken from http://stackoverflow.com/questions/3953208/value-change-listener-to-jtextfield
   * @author Boann
   */
  public static void addChangeListener(final JTextComponent text,
      final ChangeListener changeListener) {
    if (text == null || changeListener == null) {
      throw new IllegalArgumentException(
          "Both the text component and the change listener need to be non-null");
    }

    final DocumentListener dl = new DocumentListener() {
      private int lastChange = 0, lastNotifiedChange = 0;

      @Override
      public void insertUpdate(DocumentEvent e) {
        changedUpdate(e);
      }

      @Override
      public void removeUpdate(DocumentEvent e) {
        changedUpdate(e);
      }

      @Override
      public void changedUpdate(DocumentEvent e) {
        lastChange++;

        SwingUtilities.invokeLater(() -> {
          if (lastNotifiedChange != lastChange) {
            lastNotifiedChange = lastChange;
            changeListener.stateChanged(new ChangeEvent(text));
          }
        });
      }
    };

    PropertyChangeListener pcl = e -> {
      Document d1 = (Document) e.getOldValue();
      Document d2 = (Document) e.getNewValue();
      if (d1 != null) {
        d1.removeDocumentListener(dl);
      }
      if (d2 != null) {
        d2.addDocumentListener(dl);
      }
      dl.changedUpdate(null);
    };
    text.addPropertyChangeListener("document", pcl);

    Document d = text.getDocument();
    if (d != null) {
      d.addDocumentListener(dl);
    }
  }

  public static void enableComponents(Container container, boolean enabled) {
    enableComponents(container, enabled, false);
  }

  public static void enableComponents(Container container, boolean enabled,
    boolean applyToContainer) {
    enableComponents(container, enabled, applyToContainer, Collections.emptyList());
  }

  public static void enableComponents(Container container, boolean enabled,
      boolean applyToContainer, List<Component> exclusions) {
    if (applyToContainer)
      container.setEnabled(enabled);
    Component[] components = container.getComponents();
    for (Component component : components) {
      if (exclusions.contains(component)) {
        continue; // skipping excluded components
      }
      component.setEnabled(enabled);
//            if (component instanceof JScrollPane) {
//                JScrollPane jsp = (JScrollPane)component;
//                enableComponents(jsp.getViewport(), enable);
//            }
      if (component instanceof Container) {
        enableComponents((Container) component, enabled, applyToContainer);
      }
    }
  }

  /**
   * Creates a non-editable JEditorPane that has the same styling as default JLabels. Hyperlink
   * clicks are opened using the default browser.
   *
   * @param text Your text to be displayed in HTML context. Don't add the opening and closing HTML
   * tags. To include links use the regular A tags.
   */
  public static JEditorPane createClickableHtml(String text) {
    return createClickableHtml(text, true, true, null);
  }

  /**
   * Creates a non-editable JEditorPane that has the same styling as default JLabels. Hyperlink
   * clicks are opened using the default browser.
   *
   * @param applyMakeHtml Apply {@link #makeHtml(String)} function before creating Editor Pane.
   * @param text Your text to be displayed in HTML context. Don't add the opening and closing HTML
   * tags. To include links use the regular A tags.
   */
  public static JEditorPane createClickableHtml(boolean applyMakeHtml, String text) {
    return createClickableHtml(applyMakeHtml ? makeHtml(text) : text, true, true, null);
  }

  public static JScrollPane createClickableHtmlInScroll(boolean applyMakeHtml, String text) {
    return createClickableHtmlInScroll(applyMakeHtml, text, null);
  }

  public static JScrollPane createClickableHtmlInScroll(boolean applyMakeHtml, String text, Dimension preferredEditorPaneSize) {
    JEditorPane ep = createClickableHtml(applyMakeHtml ? makeHtml(text) : text, true,
        true, null);
    ep.setPreferredSize(preferredEditorPaneSize);
    JScrollPane s = new JScrollPane();
    s.setViewportView(ep);
    return s;
  }


  public static String createCssStyle() {
    return createCssStyle(null);
  }

  public static String createCssStyle(Font font) {
    // for copying style
    if (font == null) {
      font = new JLabel().getFont();
    }

    // create some css from the label's font
    StringBuilder style = new StringBuilder("font-family:" + font.getFamily() + ";");
    style.append("font-weight:").append(font.isBold() ? "bold" : "normal").append(";");
    style.append("font-size:").append(font.getSize()).append("pt;");
    return style.toString();
  }

  public static String wrapInStyledHtml(String text) {
    return wrapInStyledHtml(text, null);
  }

  public static String wrapInStyledHtml(String text, Font font) {
//    org.jsoup.nodes.Document html;
//    if (text.contains("<html")) {
//      html = Jsoup.parse(text);
//      html.body().attr("style", createCssStyle(font));
//      return html.toString();
//    } else {
//      html = org.jsoup.nodes.Document.createShell("");
//      html.body().attr("style", createCssStyle(font));
//      html.body().text(text);
//      return html.toString();
//    }
    if (text.contains("<html")) {
      text = Jsoup.parse(text).body().text();
    }

    StringBuilder sb = new StringBuilder();
    sb.append("<html><body style=\"").append(createCssStyle(font)).append("\">")
        .append(text)
        .append("</body></html>");
    return sb.toString();
  }

  public static void setJEditorPaneContent(JEditorPane ep, String text) {
    if (!"text/html".equalsIgnoreCase(ep.getContentType())) {
      ep.setText(text); // it's not styled with css in html
    } else {
      ep.setText(wrapInStyledHtml(text));
    }
  }

  public static void setJEditorPaneContent(JEditorPane ep, boolean applyMakeHtml, String text) {
    if (applyMakeHtml) {
      text = makeHtml(text);
    }
    if (!"text/html".equalsIgnoreCase(ep.getContentType())) {
      ep.setText(text); // it's not styled with css in html
    } else {
      ep.setText(wrapInStyledHtml(text));
    }
  }

//  /**
//   * Creates a non-editable JEditorPane that has the same styling as default JLabels and with
//   * hyperlinks clickable. They will be opened in the system default browser.
// * @param text Your text to be displayed in HTML context. Don't add the opening and closing HTML
//   * tags. To include links use the regular A tags.
//   * @param handleHyperlinks Add a handler for hyperlinks to be opened in the
// * @param useJlabelBackground Use default background of JLabels.
//   */
//  public static JEditorPane createClickableHtml(String text, boolean handleHyperlinks,
//      boolean useJlabelBackground) {
//    return createClickableHtml(text, handleHyperlinks, useJlabelBackground, null);
//  }

  /**
   * Creates a non-editable JEditorPane that has the same styling as default JLabels and with
   * hyperlinks clickable. They will be opened in the system default browser.
 * @param text Your text to be displayed in HTML context. Don't add the opening and closing HTML
   * tags. To include links use the regular A tags.
   * @param bgColor if {@code useJlabelBackground} is false, force this color. Can be null
   */
  public static JEditorPane createClickableHtml(String text,
      Color bgColor) {
    return createClickableHtml(text, true, false, bgColor);
  }

  /**
   * Creates a non-editable JEditorPane that has the same styling as default JLabels and with
   * hyperlinks clickable. They will be opened in the system default browser.
 * @param text Your text to be displayed in HTML context. Don't add the opening and closing HTML
   * tags. To include links use the regular A tags.
   * @param handleHyperlinks Add a handler for hyperlinks to be opened in the
 * @param useJlabelBackground Use default background of JLabels.
 * @param bgColor if {@code useJlabelBackground} is false, force this color. Can be null
   */
  public static JEditorPane createClickableHtml(String text, boolean handleHyperlinks,
      boolean useJlabelBackground, Color bgColor) {
    return createClickableHtml(text, handleHyperlinks, useJlabelBackground, bgColor,
        false);
  }

  /**
   * Creates a non-editable JEditorPane that has the same styling as default JLabels and with
   * hyperlinks clickable. They will be opened in the system default browser.
   *  @param text Your text to be displayed in HTML context. Don't add the opening and closing HTML
   * tags. To include links use the regular A tags.
   * @param handleHyperlinks Add a handler for hyperlinks to be opened in the
   * @param useJlabelBackground Use default background of JLabels.
   * @param bgColor if {@code useJlabelBackground} is false, force this color. Can be null
   * @param editable If the editor pane should be editable.
   */
  public static JEditorPane createClickableHtml(String text, boolean handleHyperlinks,
      boolean useJlabelBackground, Color bgColor, boolean editable) {

    String html1 = wrapInStyledHtml(text);
    StringBuilder sb = new StringBuilder();
    String html = sb.append("<html><body style=\"").append(createCssStyle()).append("\">")
        .append(text)
        .append("</body></html>").toString();
    if (!html.equals(html1)) {
      int a = 1;
    }
    JEditorPane ep = new JEditorPane("text/html", html1);
    ep.setEditable(editable);

    // handle link events
    if (handleHyperlinks) {
      ep.addHyperlinkListener(e -> {
        if (e.getEventType().equals(HyperlinkEvent.EventType.ACTIVATED)) {
          try {
            openBrowserOrThrow(e.getURL().toURI());
          } catch (URISyntaxException ex) {
            throw new IllegalStateException("Incorrect url/uri", ex);
          }

        }
      });
    }

    if (useJlabelBackground) {
      ep.setBackground(new JLabel().getBackground());
    } else if (bgColor != null) {
      ep.setBackground(bgColor);
    }

    return ep;
  }

  /**
   * Make the parent JDialog of a component resizable using the HierarchyListener.
   * Taken from: https://stackoverflow.com/a/7989417/88814
   */
  public static void makeDialogResizable(Component c) {
    c.addHierarchyListener(e -> {
      Window window = SwingUtilities.getWindowAncestor(c);
      if (window instanceof Dialog) {
        Dialog dialog = (Dialog)window;
        if (!dialog.isResizable()) {
          dialog.setResizable(true);
        }
      }
    });
  }

  /**
   * Tries to open the default browser.
   * @throws IllegalStateException if the operation fails.
   */
  public static void openBrowserOrThrow(String url) {
    try {
      Desktop.getDesktop().browse(new URI(url));
    } catch (IOException | URISyntaxException ex) {
      throw new IllegalStateException("Could not open link in default system browser", ex);
    }
  }


  /**
   * Tries to open the default browser.
   * @throws IllegalStateException if the operation fails.
   */
  public static void openBrowserOrThrow(URI uri) {
    try {
      Desktop.getDesktop().browse(uri);
    } catch (IOException ex) {
      throw new IllegalStateException("Could not open link in default system browser", ex);
    }
  }

  /**
   * Tries to open the default browser. Does nothing if the operation fails.
   * @param doLog Log the error with slf4j or not.
   */
  public static void openBrowserOrLog(URI uri, boolean doLog) {
    try {
      Desktop.getDesktop().browse(uri);
    } catch (IOException e) {
      if (doLog) {
        log.error("Could not open link in default system browser", e);
      }
    }
  }

  public static boolean isEnabledAndChecked(JCheckBox checkbox) {
    return checkbox.isEnabled() && checkbox.isSelected();
  }

  public static boolean isEnabledAndChecked(JToggleButton toggle) {
    return toggle.isEnabled() && toggle.isSelected();
  }

  public static Map<String, String> valuesToMap(Container origin) {
    return valuesToMap(origin, null);
  }

//  public static List<Component> allDescendants(Component comp, boolean includeOrigin) {
//    List<Component> childs = new ArrayList<>();
//    if (includeOrigin)
//    synchronized (comp.getTreeLock()) {
//      if (comp instanceof Container) {
//        Container c = (Container)comp;
//        c.
//      } else {
//        return Stream.of(comp);
//      }
//    }
//  }

  public static void traverse(Component origin, boolean includeOrigin, Consumer<Component> callback) {
    synchronized (origin.getTreeLock()) {
      ArrayDeque<Component> fifo = new ArrayDeque<>();
      if (includeOrigin)
        fifo.addLast(origin);
      else
      if (origin instanceof Container) {
        for (Component child : ((Container) origin).getComponents())
          fifo.addLast(child);
      }
      while (!fifo.isEmpty()) {
        Component comp = fifo.removeLast();
        callback.accept(comp);
        if (comp instanceof Container) {
          for (Component child : ((Container) comp).getComponents())
            fifo.addLast(child);
        }
      }
    }
  }

  /**
   * Drills down a {@link Container}, mapping all components that 1) have their name set, 2) are
   * {@link StringRepresentable} and returns the mapping.<br/>
   * Useful for persisting values from Swing windows.
   * @param compNameFilter Can be null, will accept all Component names then.
   */
  public static Map<String, String> valuesToMap(Container origin, Predicate<String> compNameFilter) {
    Map<String, Component> comps = SwingUtils.mapComponentsByName(origin, true);
    Map<String, String> map = new HashMap<>(comps.size());
    compNameFilter = compNameFilter == null ? s -> true : compNameFilter;
    for (Entry<String, Component> e : comps.entrySet()) {
      final String name = e.getKey();
      if (name == null || name.isEmpty()) {
        continue;
      }
      if (!compNameFilter.test(name)) {
        log.debug("Skipping serializing component, name filtered out: {}", name);
        continue;
      }

      final Component comp = e.getValue();
      String value;
      if (comp instanceof StringRepresentable) {
         value = ((StringRepresentable) comp).asString();
      } else if (comp instanceof JCheckBox) {
        value = Boolean.toString(((JCheckBox)comp).isSelected());
      } else if (comp instanceof JTextComponent) {
        value = ((JTextComponent)comp).getText();
      } else {
        log.debug(String
            .format("SwingUtils.valuesToMap() found component of type [%s] by name [%s] which "
                    + "does not implement [%s] and is not [%s, %s]",
                comp.getClass().getSimpleName(), comp.getName(),
                StringRepresentable.class.getSimpleName(), JCheckBox
                    .class.getSimpleName(), JTextComponent.class.getSimpleName()));
        continue;
      }

      if (value != null) {
        if (comp instanceof GhostedTextComponent && value.equals(((GhostedTextComponent) comp).getGhostText())) {
          log.debug("Skipping serializing ghost text component to map: '{}' has ghost value: '{}'", name, value);
        } else {
          map.put(name, value);
        }
      }
    }
    return map;
  }

  /**
   * Sets values for components in a {@link Container}. Components must 1) have their name set,
   * 2) be either {@link StringRepresentable} or 3) {@link JCheckBox}, {@link JTextComponent}.
   */
  public static void valuesFromMap(Container origin, Map<String, String> map) {
    Map<String, Component> comps = SwingUtils.mapComponentsByName(origin, true);
    for (Entry<String, String> kv : map.entrySet()) {
      final String name = kv.getKey();
      Component comp = comps.get(name);
      if (comp != null) {
        String s = kv.getValue();
        if (comp instanceof StringRepresentable) {
          ((StringRepresentable) comp).fromString(s);
        }
//        else if (comp instanceof JEditorPane) {
//          JEditorPane ep = (JEditorPane)comp;
//          if ("text/html".equals(ep.getContentType())) {
//            ep.setText(SwingUtils.wrapInStyledHtml(s));
//          }
//        }
        else if (comp instanceof JCheckBox) {
          ((JCheckBox)comp).setSelected(Boolean.parseBoolean(s));
        } else if (comp instanceof JTextComponent) {
          ((JTextComponent)comp).setText(s);
        } else {
          log.debug(String
              .format("SwingUtils.valuesFromMap() found component of type [%s] by name [%s] which "
                      + "does not implement [%s] and is not [%s, %s]",
                  comp.getClass().getSimpleName(), comp.getName(),
                  StringRepresentable.class.getSimpleName(), JCheckBox
                      .class.getSimpleName(), JTextComponent.class.getSimpleName()));
          continue;
        }
      }
    }
  }

  /**
   * @return Null if none of the fonts are available, otherwise the first available font found.
   */
  public static String checkFontAvailable(String... fontName) {
    final String[] fontsLocal = getAvailableFontNames();
    Font[] availableFonts = getAvailableFonts();
    for (String fontSearched : fontName) {
      for (String fontSystem : fontsLocal) {
        if (fontSystem.equals(fontSearched)) {
          return fontSystem;
        }
      }
    }
    return null;
  }

  public static String[] getAvailableFontNames() {
    String[] fontsLocal = fontNames;
    if (fontsLocal == null) {
      synchronized (fontLock) {
        fontsLocal = fontNames;
        if (fontsLocal == null) {
          GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
          fontNames = fontsLocal = ge.getAvailableFontFamilyNames();

        }
      }
    }
    return fontsLocal;
  }

  public static Font[] getAvailableFonts() {
    Font[] fontsLocal = fonts;
    if (fontsLocal == null) {
      synchronized (fontLock) {
        fontsLocal = fonts;
        if (fontsLocal == null) {
          GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
          fonts = fontsLocal = ge.getAllFonts();
        }
      }
    }
    return fontsLocal;
  }

  /**
   * Traverses from origin down the hierarchy putting all components with names set to
   * a map.
   */
  public static Map<String, Component> mapComponentsByName(Container origin,
      boolean includeOrigin) {
    if (origin == null) {
      return Collections.emptyMap();
    }
    Map<String, Component> map = new HashMap<>();
    ArrayDeque<Component> fifo = new ArrayDeque<>();
    synchronized (origin.getTreeLock()) {
      if (includeOrigin) {
        fifo.addLast(origin);
      } else {
        for (Component c : origin.getComponents()) {
          fifo.addLast(c);
        }
      }
      while (!fifo.isEmpty()) {
        Component c = fifo.removeFirst();
        String name = c.getName();
        if (!StringUtils.isNullOrWhitespace(name)) {
          map.put(name, c);
        }
        if (c instanceof Container) {
          for (Component child: ((Container)c).getComponents()) {
            fifo.addLast(child);
          }
        }
      }
    }
    return map;
  }

  /**
   * Show a message dialog wrapped into a scroll pane.
   * @param parent The parent for the dialog, null is ok.
   * @param component The component to be used as the message.
   */
  public static void showDialog(Component parent, final Component component) {
    makeDialogResizable(component);
    JOptionPane.showMessageDialog(parent, wrapInScrollForDialog(component));
  }

  /**
   * Show a message dialog wrapped into a scroll pane.
   * @param parent The parent for the dialog, null is ok.
   * @param component The component to be used as the message.
   */
  public static void showDialog(Component parent, final Component component, String title, int msgType) {
    makeDialogResizable(component);
    JOptionPane.showMessageDialog(parent, wrapInScrollForDialog(component), title, msgType);
  }

  /**
   * Show a message dialog wrapped into a scroll pane.
   * @param parent The parent for the dialog, null is ok.
   * @param component The component to be used as the message.
   */
  public static int showConfirmDialog(Component parent, final Component component) {
    makeDialogResizable(component);
    return JOptionPane.showConfirmDialog(parent, wrapInScrollForDialog(component));
  }

  public static int showChoiceDialog(Component parent, Object message, String[] options, int startingOption) {
    return JOptionPane
        .showOptionDialog(parent, message, "Delete the files?",
            JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[startingOption]);
  }

  /**
   * Wraps the given component in a scroll pane and attaches a hierarchy listener
   * that makes the parent dialog resizeable if the component is attached to a {@link Dialog}.
   * This is mainly for use with {@link JOptionPane#showMessageDialog(Component, Object)} and
   * the likes.
   */
  public static JScrollPane wrapInScrollForDialog(Component component) {
    // wrap a scrollpane around the component
    final JScrollPane scrollPane = new JScrollPane(component);
    // make the dialog resizable
    component.addHierarchyListener(e -> {
      Window window = SwingUtilities.getWindowAncestor(component);
      if (window instanceof java.awt.Dialog) {
        Dialog dialog = (Dialog) window;
        if (!dialog.isResizable()) {
          dialog.setResizable(true);
        }
      }
    });
    scrollPane.setBorder(new EmptyBorder(10, 10, 10, 10));
    return scrollPane;
  }

  /**
   * Sets the uncaught exception handler for the thread this method is invoked in
   * to a handler that shows a Swing GUI message dialog with error stacktrace.
   */
  public static void setUncaughtExceptionHandlerMessageDialog(Component parent) {
    Thread.setDefaultUncaughtExceptionHandler((t, e) -> {
      StringWriter sw = new StringWriter();
      e.printStackTrace(new PrintWriter(sw, true));
      String notes = sw.toString();

      JPanel panel = new JPanel();
      panel.setLayout(new BorderLayout());
      panel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
      panel.add(new JLabel("Something unexpected happened"), BorderLayout.PAGE_START);
      JTextArea notesArea = new JTextArea(40, 80);
      notesArea.setText(notes);
      JScrollPane notesScroller = new JScrollPane();
      notesScroller.setBorder(BorderFactory.createTitledBorder("Details: "));
      notesScroller.setViewportView(notesArea);
      panel.add(notesScroller, BorderLayout.CENTER);

      //JOptionPane.showMessageDialog(frame, "Some error details:\n\n" + notes, "Error", JOptionPane.ERROR_MESSAGE);
      //JOptionPane.showMessageDialog(frame, panel, "Error", JOptionPane.ERROR_MESSAGE);
      makeDialogResizable(panel);
      showDialog(parent, panel);
    });
  }

  /**
   * Prints the contents of the stacktrace to a string.
   */
  public static String stacktraceToString(Throwable t) {
    StringWriter sw = new StringWriter();
    t.printStackTrace(new PrintWriter(sw, true));
    return sw.toString();
  }

  /**
   * @param parent Can be null.
   */
  public static void showErrorDialogWithStacktrace(Throwable e, Component parent) {
    showErrorDialogWithStacktrace(e, parent, true);
  }

  public static String makeHtml(String html) {
    if (html == null)
      return null;
    Pattern re = Pattern.compile("^\\s*<\\s*html\\s*>\\s*");
    Pattern reNewline = Pattern.compile("(?<!<br/>)(\n)");
    Matcher m = reNewline.matcher(html);
    String s = m.replaceAll("<br/>\n");
    return re.matcher(s).find() ? s : "<html>" + s;
  }

  public static void showInfoDialog(Component parent, String htmlMessage, String title) {
    showDialog(parent, htmlMessage, title, JOptionPane.INFORMATION_MESSAGE);
  }

  public static void showWarningDialog(Component parent, String htmlMessage, String title) {
    showDialog(parent, htmlMessage, title, JOptionPane.WARNING_MESSAGE);
  }

  public static void showErrorDialog(Component parent, String htmlMessage, String title) {
    showDialog(parent, htmlMessage, title, JOptionPane.ERROR_MESSAGE);
  }

  public static void showDialog(Component parent, String htmlMessage, String title,
      int messageType) {
    JOptionPane.showMessageDialog(parent, new JLabel(makeHtml(htmlMessage)), title, messageType);
  }

  /**
   * @param parent Can be null.
   */
  public static void showErrorDialogWithStacktrace(Throwable e, Component parent, JComponent content, boolean doShowStacktrace) {
    JPanel panel = new JPanel();
    panel.setLayout(new BorderLayout());
    panel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
    if (content == null) {
      content = new JLabel("Something unexpected happened (" + e.getClass().getSimpleName() + ")");
    }
    panel.add(content, BorderLayout.PAGE_START);
    JTextArea notesArea = new JTextArea(40, 80);
    if (doShowStacktrace) {
      notesArea.setText(stacktraceToString(e));
    } else {
      notesArea.setText(e.getMessage());
    }
    JScrollPane notesScroller = new JScrollPane();
    notesScroller.setBorder(BorderFactory.createTitledBorder("Details: "));
    notesScroller.setViewportView(notesArea);
    panel.add(notesScroller, BorderLayout.CENTER);

    //JOptionPane.showMessageDialog(frame, "Some error details:\n\n" + notes, "Error", JOptionPane.ERROR_MESSAGE);
    //JOptionPane.showMessageDialog(frame, panel, "Error", JOptionPane.ERROR_MESSAGE);
    makeDialogResizable(panel);
    showDialog(parent, panel);
  }


  /**
   * @param parent Can be null.
   */
  public static void showErrorDialogWithStacktrace(Throwable e, Component parent, boolean doShowStacktrace) {
    showErrorDialogWithStacktrace(e, parent, null, doShowStacktrace);
  }

  /**
   * @param path If the passed path already exists, just returns it.
   * @return null if no existing Path could be found on the filesystem all the way up to root.
   */
  public static Path findExistingUpstreamPath(Path path) {
    if (path == null || Files.exists(path)) {
      return path;
    } else {
      return findExistingUpstreamPath(path.getParent());
    }
  }

  public static JFrame findParentFrame(Component origin) {
    Component parentFrameForDialog = findParentFrameForDialog(origin);
    if (parentFrameForDialog instanceof JFrame) {
      return (JFrame) parentFrameForDialog;
    }
    return null;
  }

  /**
   * Bubbles up the component hierarchy searching for first instance of a {@link JFrame}.
   */
  public static Component findParentFrameForDialog(Component origin) {
    if (origin == null) {
      return null;
    }
    if (origin instanceof JFrame) {
      return origin;
    }

    Container parent = origin.getParent();
    while (parent != null && !(parent instanceof JFrame)) {
      parent = parent.getParent();
    }
    return parent;
  }

  /**
   * Tries to set the LAF to native for the platform. Does nothing if the LAF is not available.
   *
   * @return true if setting the LAF succeeded.
   */
  public static boolean setPlatformLookAndFeel() {
    try {
      String laf = UIManager.getSystemLookAndFeelClassName();
      UIManager.setLookAndFeel(laf);
      return true;
    } catch (Exception e) {
      return false;
    }
  }

  /**
   * Tries to set the LAF to native for the platform or Nimbus if failed.
   *
   * @return true if setting the LAF succeeded.
   */
  public static boolean setPlatformLafOrNimbus() {
    if (setPlatformLookAndFeel())
      return true;
    try {
      for (UIManager.LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {
        if ("Nimbus".equals(info.getName())) {
          UIManager.setLookAndFeel(info.getClassName());
          return true;
        }
      }
    } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException ignored) {}
    return false;
  }

  /**
   * Centers a JFrame on screen.
   *
   * @param frame the frame to be centered
   */
  public static void centerFrame(JFrame frame) {
    Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
    frame.setLocation(dim.width / 2 - frame.getSize().width / 2,
        dim.height / 2 - frame.getSize().height / 2);
  }

  /**
   * Sets the icons for a frame. These are also used to display icons in the taskbar.
   *
   * @param frame The frame to set the icons for.
   * @param iconPaths The simplest way is to provide just the file names.
   * @param classToFindResources A class relative to which the icons will be searched. This is a
   * kludge to make things more fool-proof.
   */
  public static void setFrameIcons(JFrame frame, java.util.List<String> iconPaths,
      Class<?> classToFindResources) {
    java.util.List<Image> icons = new ArrayList<>();
    for (String iconPath : iconPaths) {
      java.net.URL imgURL = classToFindResources.getResource(iconPath);
      ImageIcon image = new ImageIcon(imgURL);
      icons.add(image.getImage());
    }
    frame.setIconImages(icons);
  }

  /**
   * Nicely closes the frame, respecting its {@code setOnCloseOperation()} settings.
   *
   * @param frame the frame to be closed
   */
  public static void closeFrameNicely(JFrame frame) {
    frame.dispatchEvent(new WindowEvent(frame, WindowEvent.WINDOW_CLOSING));
  }

  public static boolean isGraphicalEnvironmentAvailable() {
    boolean headless = true;

    String nm = System.getProperty("java.awt.headless");

    if (nm == null) {
      /* No need to ask for DISPLAY when run in a browser */
      if (System.getProperty("javaplugin.version") != null) {
        headless = Boolean.FALSE;
      } else {
        String osName = System.getProperty("os.name");
        headless = ("Linux".equals(osName) ||
            "SunOS".equals(osName)) && (System.getenv("DISPLAY") == null);
      }
    } else if (nm.equals("true")) {
      headless = Boolean.TRUE;
    } else {
      headless = Boolean.FALSE;
    }

    return !headless;
  }

  public static void userShowDialog(Component frame, final Component component) {
    // wrap a scrollpane around the component
    JScrollPane scrollPane = new JScrollPane(component);
    // make the dialog resizable
    component.addHierarchyListener(e -> {
      Window window = SwingUtilities.getWindowAncestor(component);
      if (window instanceof Dialog) {
        Dialog dialog = (Dialog) window;
        if (!dialog.isResizable()) {
          dialog.setResizable(true);
        }
      }
    });
    // display them in a message dialog
    JOptionPane.showMessageDialog(frame, scrollPane);
  }

  public static void userShowError(Component frame, String stacktrace) {
    JPanel panel = new JPanel();
    panel.setLayout(new BorderLayout());
    panel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
    panel.add(new JLabel("Something unexpected happened"), BorderLayout.PAGE_START);
    JTextArea notesArea = new JTextArea(40, 80);
    notesArea.setText(stacktrace);
    JScrollPane notesScroller = new JScrollPane();
    notesScroller.setBorder(BorderFactory.createTitledBorder("Details: "));
    notesScroller.setViewportView(notesArea);
    panel.add(notesScroller, BorderLayout.CENTER);
    //JOptionPane.showMessageDialog(frame, "Some error details:\n\n" + notes, "Error", JOptionPane.ERROR_MESSAGE);
    //JOptionPane.showMessageDialog(frame, panel, "Error", JOptionPane.ERROR_MESSAGE);
    makeDialogResizable(panel);
    userShowDialog(frame, panel);
  }
}
