//
// Created by James Brandenburg on 6/4/15.
//

#ifndef ROOBARB_XMLCANVAS_H
#define ROOBARB_XMLCANVAS_H

// STL
#include <map>
using namespace std;

#include "Logger.h"
#include "XmlConfig.h"
#include "RooBarb/XmlPad.h"

#include "TCanvas.h"


namespace jdb{
    class XmlCanvas : public IObject{

    protected:
        string name, title;

        int pxWidth, pxHeight;
        int nCol, nRow;

        TCanvas * rootCanvas = nullptr;
        map<string, shared_ptr<XmlPad> > pads;

    public:
    virtual const char * classname() const { return "XmlCanvas"; }
        string getName() { return name; }

        XmlCanvas( XmlConfig &cfg, string _nodePath ) {
            if ( cfg.exists( _nodePath ) ) {
                string preNode = cfg.cn(_nodePath);

                pxWidth = cfg.getInt( ":width", cfg.getInt( ":w", 800 ) );
                pxHeight = cfg.getInt( ":height", cfg.getInt( ":h", 800 ) );

                name = cfg.getString( ":name", cfg.getString( ":n", "c" ) );
                title = cfg.getString( ":title", cfg.getString( ":t", "c" ) );

                nCol = cfg.getInt( ":columns", cfg.getInt( ":cols", cfg.getInt( ":nCols", 12 ) ) );
                nRow = cfg.getInt( ":rows", 12 );

                DEBUGC( "name=" << name << ", title=" << title <<", width=" << pxWidth << ", height=" << pxHeight );
                rootCanvas = new TCanvas( name.c_str(), title.c_str(), pxWidth, pxHeight );
                cfg.cn(preNode);
                
                createPads( cfg, _nodePath );

                
            } else {
                TRACEC("No XmlConfig Given : Creating Default Canvas" );
                rootCanvas = new TCanvas( "XmlCanvas", "XmlCanvas", 800, 1200 );
            }

        }

        void cd( int iPad = 0 ) {
            if ( !rootCanvas )
                return;
            rootCanvas->cd(iPad);
        }

        void saveImage( string iname ) {
            if ( !rootCanvas )
                return;

            rootCanvas->Print( iname.c_str() );
        }

        shared_ptr<XmlPad> activatePad( string padName ){
            LOG_F( INFO, "activate pad: %s", padName.c_str()  );
            
            if ( pads.count( padName ) > 0 ){
                LOG_F( INFO, "Found: %s", padName.c_str()  );
                shared_ptr<XmlPad> xpad =  pads[ padName ];
                if ( nullptr != xpad ){
                    xpad->cd();
                } else {
                    LOG_F( INFO, "PAD IS NULL" );
                }
                return xpad;
            }
            
            return nullptr;
        }

        TCanvas * getCanvas(){
            return rootCanvas;
        }

        

    protected:
        void createPads( XmlConfig &cfg, string _nodePath ){
            LOG_F( INFO, "Creating pads" );
            vector<string> children = cfg.childrenOf( _nodePath, "Pad" );
            LOG_S(INFO) << "Found " << children.size();
            for ( string path : children ){
                LOG_S(INFO) << "Creating Pad From " << path;
                this->cd();
                // cfg.cn( path );
                string name = cfg.getString( path + ":name", cfg.getString( path+":n", "" ) );
                LOG_S(INFO) << "Creating Pad named " << name;
                if ( "" != name ){
                    pads[ name ] = shared_ptr<XmlPad>( new XmlPad( cfg, path, nRow, nCol ) );

                }
                
            }

        }

    };
}// namespace

#endif //ROOBARB_XMLCANVAS_H
